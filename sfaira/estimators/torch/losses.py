import numpy as np
import torch


class LossLoglikelihoodNb(torch.nn.Module):

    def __init__(self, average=True):
        super(LossLoglikelihoodNb, self).__init__()
        self.average = average

    def forward(self, preds, target):
        """Implements the negative log likelihood loss as VAE reconstruction loss"""
        x = target
        loc, scale = torch.chunk(preds, chunks=2, dim=1)

        eta_loc = torch.log(loc)
        eta_scale = torch.log(scale)

        log_r_plus_mu = torch.log(scale + loc)

        ll = torch.lgamma(scale + x)
        ll = ll - torch.lgamma(x + torch.ones_like(x))
        ll = ll - torch.lgamma(scale)
        ll = ll + torch.multiply(x, eta_loc - log_r_plus_mu) + torch.multiply(scale, eta_scale - log_r_plus_mu)

        ll = torch.clamp(ll, min=-300, max=300)
        neg_ll = -ll
        if self.average:
            neg_ll = torch.mean(neg_ll)
        else:
            # sum over features, average over batch
            neg_ll = neg_ll.sum(dim=1).sum(dim=1)
        return neg_ll


class LossLoglikelihoodGaussian(torch.nn.Module):

    def __init__(self, average=True):
        super(LossLoglikelihoodGaussian, self).__init__()
        self.average = average

    def forward(self, preds, target):
        """Implements the gaussian log likelihood loss as VAE reconstruction loss"""
        loc, scale = torch.chunk(preds, chunks=2, dim=1)

        ll = -torch.log(scale * torch.sqrt(2. * np.pi)) - 0.5 * torch.square((target - loc) / scale)
        ll = torch.clamp(ll, min=-300, max=300)
        neg_ll = -ll
        if self.average:
            neg_ll = torch.mean(neg_ll)
        else:
            # sum over features, average over batch
            neg_ll = torch.mean(torch.sum(neg_ll, dim=1), dim=0)
        return neg_ll


class LossCrossentropyAgg(torch.nn.Module):

    def __init__(self):
        super(LossCrossentropyAgg, self).__init__()

    def forward(self, preds, target):
        """ Modified crossentropy that aggregates allowed output classes into single class. """
        preds = torch.clamp(preds, min=1e-10, max=1.)
        ll_cce_agg = -torch.log(torch.mean(target * preds, dim=1, keepdim=False))
        return ll_cce_agg


class KLLoss(torch.nn.Module):

    def __init__(self):
        super(KLLoss, self).__init__()
        self.beta = self.register_buffer('beta', torch.Tensor(1.))

    def forward(self, preds, target):
        expected_logqz_x, expected_logpz = torch.chunk(preds, chunks=2, dim=1)

        kl_loss = torch.mean(expected_logqz_x - expected_logpz, dim=0)
        return self.beta * kl_loss
