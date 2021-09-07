import numpy as np
try:
    import tensorflow as tf
except ImportError:
    tf = None
from typing import Union, Tuple

from sfaira.models.embedding.output_layers import NegBinOutput, NegBinSharedDispOutput, NegBinConstDispOutput, \
    GaussianOutput, GaussianSharedStdOutput, GaussianConstStdOutput
from sfaira.versions.topologies import TopologyContainer
from sfaira.models.embedding.base import BasicModelKerasEmbedding
from sfaira.models.pp_layer import PreprocInput
from sfaira.models.made import MaskingDense


class Sampling(tf.keras.layers.Layer):
    """Uses (z_mean, z_log_var) to sample z."""

    def call(self, inputs, **kwargs):
        z_mean, z_log_var = inputs
        batch = tf.shape(z_mean)[0]
        dim = tf.shape(z_mean)[1]
        epsilon = tf.keras.backend.random_normal(shape=(batch, dim))
        return z_mean + tf.exp(0.5 * z_log_var) * epsilon


class Encoder(tf.keras.layers.Layer):
    """Maps input to embedding space"""

    def __init__(
            self,
            latent_dim: Tuple,
            dropout_rate,
            l1_coef: float,
            l2_coef: float,
            batchnorm: bool,
            activation: str,
            kernel_initializer: str,
            name='encoder',
            **kwargs
    ):
        super().__init__(name=name, **kwargs)

        self.fwd_pass = []
        for i, hid_size in enumerate(latent_dim[:-1]):
            self.fwd_pass.append(
                tf.keras.layers.Dense(
                    units=hid_size,
                    activation=activation,
                    kernel_initializer=kernel_initializer,
                    kernel_regularizer=tf.keras.regularizers.l1_l2(l1=l1_coef, l2=l2_coef),
                    name='enc%s' % i
                )
            )
            if batchnorm:
                self.fwd_pass.append(
                    tf.keras.layers.BatchNormalization(
                        center=True,
                        scale=True
                    )
                )
            self.fwd_pass.append(
                tf.keras.layers.Dropout(
                    dropout_rate,
                    noise_shape=None,
                    seed=None
                )
            )

        # final layer
        self.dense_mean = tf.keras.layers.Dense(
            units=latent_dim[-1],
            activation="linear",
            name='z_mean'
        )
        self.dense_log_sigma = tf.keras.layers.Dense(
            units=latent_dim[-1],
            activation="linear",
            name='z_log_var'
        )
        # IAF layers
        self.dense_h = tf.keras.layers.Dense(
            units=latent_dim[-1],
            activation='sigmoid',
            name='h'
        )

    def call(self, inputs, **kwargs):
        x = inputs
        for layer in self.fwd_pass:
            x = layer(x)
        # final layer
        z_mean = self.dense_mean(x)
        z_log_var = self.dense_log_sigma(x)
        h = self.dense_h(x)

        return z_mean, z_log_var, h


class IAF(tf.keras.layers.Layer):
    def __init__(
        self,
        bottleneck: int,
        n_iaf: int,
        l1_coef: float,
        l2_coef: float,
        masking_dim=320,
        n_made=2,
        activation="relu",
        name='iaf',
        **kwargs
    ):
        """
        Transforms latent space with simple distribution to one with a more flexible one.

        :param latent_dim: dimension of bottleneck
        :param masking_dim: units of autoregressive networks
        :param n_made: number of autoregressive networks
        :param n_iaf: number of IAF transforms
        :param l1_coef:
        :param l2_coef:
        :param activation:
        :param name:
        :param kwargs:
        """
        super().__init__(name=name, **kwargs)
        self.m_fwd_pass = []
        self.s_fwd_pass = []
        for i in range(n_iaf):
            self.m_fwd_pass.append(
                MaskingDense(
                    units=masking_dim,
                    out_units=bottleneck,
                    hidden_layers=n_made,
                    activation=activation,
                    out_activation='linear',
                    # batchnorm=True,
                    name='mean_%s' % i)
            )
            self.s_fwd_pass.append(
                MaskingDense(
                    units=masking_dim,
                    out_units=bottleneck,
                    hidden_layers=n_made,
                    kernel_regularizer=tf.keras.regularizers.l1_l2(l1=l1_coef, l2=l2_coef),
                    out_bias_initializer='ones',
                    activation=activation,
                    out_activation='sigmoid',
                    # batchnorm=True,
                    name='sigma_%s' % i)
            )

    def call(self, inputs, **kwargs):
        """
        :param inputs: initial latent space z and dense layer h
        :param kwargs:
        :return: z_t: transformed latent space
        :return: s_t_sigmas: list of applied transformations
        """
        z_t, h = inputs

        s_t_sigmas = []
        for i, (m_iaf, s_iaf) in enumerate(zip(self.m_fwd_pass, self.s_fwd_pass)):
            if i > 0:
                z_t = tf.keras.layers.Lambda(lambda x: tf.keras.backend.reverse(x, axes=-1))(z_t)
                # z_t = tf.gather(z_t, tf.random.shuffle(tf.range(tf.shape(z_t)[0])))  # random permutation
            m_t = m_iaf([z_t, h])
            sigma_t = s_iaf([z_t, h])
            z_t = z_t * sigma_t + m_t * (1. - sigma_t)
            s_t_sigmas.append(sigma_t)
        return z_t, s_t_sigmas


class Decoder(tf.keras.layers.Layer):
    """Maps latent space sample back to output"""

    def __init__(
            self,
            latent_dim: Tuple,
            dropout_rate,
            l1_coef: float,
            l2_coef: float,
            batchnorm: bool,
            activation: str,
            kernel_initializer: str,
            name='decoder',
            **kwargs
    ):
        super().__init__(name=name, **kwargs)

        self.fwd_pass = []
        for i, hid_size in enumerate(latent_dim):
            self.fwd_pass.append(
                tf.keras.layers.Dense(
                    units=hid_size,
                    activation=activation,
                    kernel_initializer=kernel_initializer,
                    kernel_regularizer=tf.keras.regularizers.l1_l2(l1=l1_coef, l2=l2_coef),
                    name='dec%s' % i
                )
            )
            if batchnorm:
                self.fwd_pass.append(
                    tf.keras.layers.BatchNormalization(
                        center=True,
                        scale=True
                    )
                )
            self.fwd_pass.append(
                tf.keras.layers.Dropout(
                    dropout_rate,
                    noise_shape=None,
                    seed=None
                )
            )

    def call(self, inputs, **kwargs):
        x = inputs
        for layer in self.fwd_pass:
            x = layer(x)
        return x


class ModelKerasVaeIAF(BasicModelKerasEmbedding):

    def __init__(
            self,
            in_dim,
            latent_dim=(128, 64, 2, 64, 128),
            n_iaf=2,
            dropout_rate=0.1,
            l2_coef: float = 0.,
            l1_coef: float = 0.,
            mc_samples=10,
            batchnorm=False,
            activation='tanh',
            init='glorot_uniform',
            output_layer="nb"
    ):
        super(ModelKerasVaeIAF, self).__init__()
        # Check length of latent dim to divide encoder-decoder stack:
        if len(latent_dim) % 2 == 1:
            n_layers_enc = len(latent_dim) // 2 + 1
        else:
            raise ValueError("len(latent_dim)=%i should be uneven to provide a defined bottleneck" % len(latent_dim))

        inputs_encoder = tf.keras.Input(shape=(in_dim,), name='counts')
        inputs_sf = tf.keras.Input(shape=(1,), name='size_factors')
        inputs_encoder_pp = PreprocInput()(inputs_encoder)
        encoder = Encoder(
            latent_dim=latent_dim[:n_layers_enc],
            dropout_rate=dropout_rate,
            l1_coef=l1_coef,
            l2_coef=l2_coef,
            batchnorm=batchnorm,
            activation=activation,
            kernel_initializer=init
        )
        iaf = IAF(
            bottleneck=latent_dim[n_layers_enc - 1],
            n_iaf=n_iaf,
            l1_coef=l1_coef,
            l2_coef=l2_coef
        )
        decoder = Decoder(
            latent_dim=latent_dim[n_layers_enc:],
            dropout_rate=dropout_rate,
            l1_coef=l1_coef,
            l2_coef=l2_coef,
            batchnorm=batchnorm,
            activation=activation,
            kernel_initializer=init
        )
        sampling = Sampling()
        z_mean, z_log_var, h = encoder(inputs_encoder_pp)

        # for prediction
        z_0 = sampling([z_mean, z_log_var])
        z_t, s_t_sigmas = iaf([z_0, h])
        z_t_square_mc = 0
        z_t_mean = 0

        # for kl divergence
        for i in range(mc_samples):
            z = sampling([z_mean, z_log_var])
            z, s_t_sigmas = iaf([z, h])
            z_t_square_mc += tf.square(z)
            z_t_mean += z
        z_t_square_mc = z_t_square_mc / mc_samples
        z_t_mean = z_t_mean / mc_samples

        cum_s_t_log_var = 0
        for s_t_sigma in s_t_sigmas:
            cum_s_t_log_var += 2.0 * tf.keras.backend.log(s_t_sigma + 1e-10)

        expected_logqz_x = -0.5 * tf.reduce_sum(1 + z_log_var + cum_s_t_log_var, axis=1)
        expected_logpz = -0.5 * tf.reduce_sum(z_t_square_mc, axis=1)

        expected_densities = tf.keras.layers.Concatenate(axis=1, name='kl')([
            tf.expand_dims(expected_logqz_x, axis=1),
            tf.expand_dims(expected_logpz, axis=1)])

        output_decoder = decoder(z_t)

        if output_layer == 'nb':
            output_decoder_expfamily = NegBinOutput(original_dim=in_dim)((output_decoder, inputs_sf))
        elif output_layer == 'nb_shared_disp':
            output_decoder_expfamily = NegBinSharedDispOutput(original_dim=in_dim)((output_decoder, inputs_sf))
        elif output_layer == 'nb_const_disp':
            output_decoder_expfamily = NegBinConstDispOutput(original_dim=in_dim)((output_decoder, inputs_sf))
        elif output_layer == 'gaussian':
            output_decoder_expfamily = GaussianOutput(original_dim=in_dim)((output_decoder, inputs_sf))
        elif output_layer == 'gaussian_shared_disp':
            output_decoder_expfamily = GaussianSharedStdOutput(original_dim=in_dim)((output_decoder, inputs_sf))
        elif output_layer == 'gaussian_const_disp':
            output_decoder_expfamily = GaussianConstStdOutput(original_dim=in_dim)((output_decoder, inputs_sf))
        else:
            raise ValueError("tried to access a non-supported output layer %s" % output_layer)
        output_decoder_expfamily_concat = tf.keras.layers.Concatenate(axis=1, name="neg_ll")(output_decoder_expfamily)

        self.encoder_model = tf.keras.Model(
            inputs=[inputs_encoder, inputs_sf],
            outputs=[z_t, z_t_mean, z_0],
            name="encoder_model"
        )
        self.training_model = tf.keras.Model(
            inputs=[inputs_encoder, inputs_sf],
            outputs=[output_decoder_expfamily_concat, expected_densities],
            name="autoencoder"
        )

    def predict_embedding(self, x, variational=False, return_z0=False):
        if return_z0 and variational:
            z_t, z_t_mean, z_0 = self.encoder_model.predict(x)
            return z_t, z_t_mean, z_0
        elif not return_z0 and variational:
            z_t, z_t_mean = self.encoder_model.predict(x)[:2]
            return z_t, z_t_mean
        elif return_z0 and not variational:
            z_t_mean, z_0 = self.encoder_model.predict(x)[1:]
            return z_t_mean, z_0
        elif not return_z0 and not variational:
            z_t_mean = self.encoder_model.predict(x)[1]
            return z_t_mean


class ModelVaeIAFVersioned(ModelKerasVaeIAF):
    def __init__(
            self,
            topology_container: TopologyContainer,
            override_hyperpar: Union[dict, None] = None
    ):
        hyperpar = topology_container.topology["hyper_parameters"]
        if override_hyperpar is not None:
            for k in list(override_hyperpar.keys()):
                hyperpar[k] = override_hyperpar[k]
        super().__init__(
            in_dim=topology_container.n_var,
            **hyperpar
        )
        print('passed hyperpar: \n', hyperpar)
        self._topology_id = topology_container.topology_id
        self.genome_size = topology_container.n_var
        self.model_class = "embedding"
        self.model_type = topology_container.model_type
        self.hyperparam = dict(
            list(hyperpar.items()) +  # noqa: W504
            [
                ("topology_id", self._topology_id),
                ("genome_size", self.genome_size),
                ("model_class", self.model_class),
                ("model_type", self.model_type)
            ]
        )
