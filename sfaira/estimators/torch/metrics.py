import numpy as np
import torch
import torchmetrics


def custom_mse(y_true, y_pred, sample_weight=None):
    y_pred = torch.chunk(y_pred, chunks=2, dim=1)[0]
    se = torch.square(torch.subtract(y_true, y_pred))
    se_red = torch.mean(se)
    return se_red


def custom_mean_squared_logp1_error(y_true, y_pred, sample_weight=None):
    y_pred = torch.chunk(y_pred, chunks=2, dim=1)[0]
    y_true = torch.log(y_true + 1.)
    y_pred = torch.log(y_pred + 1.)
    se = torch.square(torch.subtract(y_true, y_pred))
    se_red = torch.mean(se)
    return se_red


def custom_negll_nb(y_true, y_pred, sample_weight=None):
    x = y_true
    loc, scale = torch.chunk(y_pred, chunks=2, dim=1)

    eta_loc = torch.log(loc)
    eta_scale = torch.log(scale)

    log_r_plus_mu = torch.log(scale + loc)

    ll = torch.lgamma(scale + x)
    ll = ll - torch.lgamma(x + torch.ones_like(x))
    ll = ll - torch.lgamma(scale)
    ll = ll + torch.multiply(x, eta_loc - log_r_plus_mu) + torch.multiply(scale, eta_scale - log_r_plus_mu)

    ll = torch.clamp(ll, min=-300, max=300)
    neg_ll = -ll
    neg_ll = torch.mean(neg_ll)
    return neg_ll


def custom_negll_gaussian(y_true, y_pred, sample_weight=None):
    loc, scale = torch.chunk(y_pred, chunks=2, dim=1)

    ll = -torch.log(scale * torch.sqrt(2. * np.pi)) - 0.5 * torch.square((y_true - loc) / scale)
    ll = torch.clamp(ll, min=-300, max=300)
    neg_ll = -ll
    neg_ll = torch.mean(neg_ll)
    return neg_ll


def custom_kl(y_true, y_pred, sample_weight=None):
    expected_logqz_x, expected_logpz = torch.chunk(y_pred, chunks=2, dim=1)
    kl_loss = torch.mean(expected_logqz_x - expected_logpz)
    return kl_loss


def custom_cce_agg(y_true, y_pred, sample_weight=None):
    y_pred = torch.clamp(y_pred, min=1e-10, max=1.)
    ll_cce_agg = -torch.log(torch.mean(y_true * y_pred, dim=1, keepdim=False))
    return ll_cce_agg


class CustomAccAgg(torchmetrics.Metric):

    def __init__(self, dist_sync_on_step=False):
        super(CustomAccAgg, self).__init__(dist_sync_on_step=dist_sync_on_step)
        self.add_state('acc', default=torch.Tensor([0.]))
        self.add_state('count', default=torch.Tensor([0.]))

    def update(self, preds, target, sample_weight=None):
        phat_pos_agg = torch.sum(target * preds, dim=1, keepdim=True)
        acc_agg = (phat_pos_agg > torch.max((torch.ones_like(target) - target) * preds, dim=1)).type(target.dtype)
        # Do not use weighting for accuracy.
        self.acc += self.acc(torch.mean(acc_agg))
        self.count += 1.

    def compute(self):
        return torch.divide(self.acc, self.count)


class CustomTprClasswise(torchmetrics.Metric):

    def __init__(self, k: int, dist_sync_on_step=False):
        super(CustomTprClasswise, self).__init__(dist_sync_on_step=dist_sync_on_step)
        self.add_state('tp', default=torch.zeros(size=(k, )))
        self.add_state('fn', default=torch.zeros(size=(k, )))

    def update(self, preds, target, sample_weight=None):
        tp_by_class = torch.sum((preds == torch.max(preds, dim=1, keepdim=True)).type(target.dtype) * target, dim=0)
        fn_by_class = torch.sum((preds < torch.max(preds, dim=1, keepdim=True)).type(target.dtype) * target, dim=0)
        self.tp += tp_by_class
        self.fn += fn_by_class

    def compute(self):
        # Catch division by zero, in that case tpr becomes zero after clipping the divisor to 1.
        divisor = torch.clamp(self.tp + self.fn, min=1., max=np.inf)
        tpr = torch.mean(self.tp / divisor)
        return tpr


class CustomFprClasswise(torchmetrics.Metric):

    def __init__(self, k: int, dist_sync_on_step=False):
        super(CustomFprClasswise, self).__init__(dist_sync_on_step=dist_sync_on_step)
        self.add_state('fp', default=torch.zeros(size=(k, )))
        self.add_state('tn', default=torch.zeros(size=(k, )))

    def update(self, preds, target, sample_weight=None):
        fp_by_class = torch.sum(
            (preds == torch.max(preds, dim=1, keepdim=True)).type(target.dtype) * (1. - target),
            dim=0
        )
        tn_by_class = torch.sum(
            (preds < torch.max(preds, dim=1, keepdim=True)).type(target.dtype) * (1. - target),
            dim=0
        )
        self.fp += fp_by_class
        self.tn += tn_by_class

    def compute(self):
        # Catch division by zero, in that case fpr becomes zero after clipping the divisor to 1.
        divisor = torch.clamp(self.fp + self.tn, min=1., max=np.inf)
        fpr = torch.mean(self.fp / divisor)
        return fpr


class CustomF1Classwise(torchmetrics.Metric):

    def __init__(self, k: int, dist_sync_on_step=False):
        super(CustomF1Classwise, self).__init__(dist_sync_on_step=dist_sync_on_step)
        self.add_state('tp', default=torch.zeros(size=(k, )))
        self.add_state('fp', default=torch.zeros(size=(k, )))
        self.add_state('fn', default=torch.zeros(size=(k, )))

    def update(self, preds, target, sample_weight=None):
        tp_by_class = torch.sum((
            preds == torch.max(preds, dim=1, keepdim=True)).type(target.dtype) * target, dim=0)
        fp_by_class = torch.sum(
            (preds == torch.max(preds, dim=1, keepdim=True)).type(target.dtype) * (1. - target),
            dim=0
        )
        fn_by_class = torch.sum((preds < torch.max(preds, dim=1, keepdim=True)).type(target.dtype) * target, dim=0)
        self.tp.assign_add(tp_by_class)
        self.fp.assign_add(fp_by_class)
        self.fn.assign_add(fn_by_class)

    def compute(self):
        # Catch divisions by zero, in that case precision or recall become zero after clipping the divisor to 1.
        divisor_precision = torch.clamp(self.tp + self.fp, min=1, max=np.inf)
        divisor_recall = torch.clamp(self.tp + self.fn, min=1., max=np.inf)
        precision = self.tp / divisor_precision
        recall = self.tp / divisor_recall
        precision = torch.clamp(precision, min=1e-100, max=np.inf)
        recall = torch.clamp(recall, min=1e-100, max=np.inf)
        f1 = torch.mean(2 * 1 / (1 / precision + 1 / recall))
        return f1
