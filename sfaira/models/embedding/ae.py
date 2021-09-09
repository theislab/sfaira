import numpy as np
try:
    import tensorflow as tf
except ImportError:
    tf = None
from typing import List, Union, Tuple

from sfaira.models.embedding.output_layers import NegBinOutput, NegBinSharedDispOutput, NegBinConstDispOutput, \
    GaussianOutput, GaussianSharedStdOutput, GaussianConstStdOutput
from sfaira.versions.topologies import TopologyContainer
from sfaira.models.embedding.base import BasicModelKerasEmbedding
from sfaira.models.pp_layer import PreprocInput


class Encoder(tf.keras.layers.Layer):
    """Project input data to latent embedding space"""

    def __init__(
            self,
            latent_dim: Tuple,
            l2_coef: float,
            l1_coef: float,
            hidden_dropout: List[float],
            input_dropout: float,
            batchnorm: bool,
            activation: str,
            kernel_initializer: str,
            name='encoder',
            **kwargs
    ):
        super().__init__(name=name, **kwargs)

        self.layer_list = []

        if input_dropout > 0.0:
            self.layer_list.append(tf.keras.layers.Dropout(input_dropout, name='input_dropout'))

        # build layer list
        for i, (hid_size, hid_drop) in enumerate(zip(latent_dim, hidden_dropout)):
            # At the bottleneck layer:
            # 1. a linear activation function is used to not restrict the support of the hidden units
            # Additionally, before and after the bottleneck dense layer:
            # 2. no batch normalisation is used to not scale the variance of inactive units up and the activity of
            # active units down.
            # 3. no dropout is used to not confound the required bottleneck dimension with dropout rate.
            self.layer_list.append(
                tf.keras.layers.Dense(
                    hid_size,
                    activation=activation if i < (len(latent_dim) - 1) else "linear",
                    kernel_initializer=kernel_initializer,
                    kernel_regularizer=tf.keras.regularizers.l1_l2(l1_coef, l2_coef),
                    name='enc%s' % i
                )
            )
            if batchnorm and i < (len(latent_dim) - 2):
                self.layer_list.append(tf.keras.layers.BatchNormalization(center=True, scale=False))
            if hid_drop > 0.0 and i < (len(latent_dim) - 2):
                if activation == "selu":
                    self.layer_list.append(tf.keras.layers.AlphaDropout(hid_drop, name='enc_%s_drop' % i))
                else:
                    self.layer_list.append(tf.keras.layers.Dropout(hid_drop, name='enc_%s_drop' % i))

    def call(self, x, **kwargs):
        for layer in self.layer_list:
            x = layer(x)
        return x


class Decoder(tf.keras.layers.Layer):
    """Recover output from latent embedding"""

    def __init__(
            self,
            latent_dim: Tuple,
            l2_coef: float,
            l1_coef: float,
            hidden_dropout: List[float],
            batchnorm: bool,
            activation: str,
            kernel_initializer: str,
            name='decoder',
            **kwargs
    ):
        super().__init__(name=name, **kwargs)

        self.layer_list = []

        # build layer list
        for i, (hid_size, hid_drop) in enumerate(zip(latent_dim, hidden_dropout)):
            self.layer_list.append(
                tf.keras.layers.Dense(
                    hid_size,
                    activation=activation,
                    kernel_initializer=kernel_initializer,
                    kernel_regularizer=tf.keras.regularizers.l1_l2(l1_coef, l2_coef),
                    name='dec%s' % i
                )
            )
            if batchnorm:
                self.layer_list.append(
                    tf.keras.layers.BatchNormalization(center=True, scale=False)
                )

            if hid_drop > 0.0:
                if activation == "selu":
                    self.layer_list.append(tf.keras.layers.AlphaDropout(hid_drop, name='dec_%s_drop' % i))
                else:
                    self.layer_list.append(tf.keras.layers.Dropout(hid_drop, name='dec_%s_drop' % i))

    def call(self, x, **kwargs):
        for layer in self.layer_list:
            x = layer(x)
        return x


class ModelKerasAe(BasicModelKerasEmbedding):
    """Combines the encoder and decoder into an end-to-end model for training."""
    # Note: Original DCA implementation uses l1_l2 regularisation also on last layer (nb) - missing here
    # Note: Original DCA implementation uses softplus function instead of exponential as dispersion activation

    def __init__(
            self,
            in_dim,
            latent_dim: Tuple = (64, 32, 64),
            l2_coef: float = 0.,
            l1_coef: float = 0.,
            dropout_rate: float = 0.,
            input_dropout: float = 0.,
            batchnorm: bool = True,
            activation='relu',
            init='glorot_uniform',
            output_layer="nb"
    ):
        # Check length of latent dim to divide encoder-decoder stack:
        if len(latent_dim) % 2 == 1:
            n_layers_enc = len(latent_dim) // 2 + 1
        else:
            raise ValueError("len(latent_dim)=%i should be uneven to provide a defined bottleneck" % len(latent_dim))

        self.in_dim = in_dim
        self.latent_dim = latent_dim
        out_dim = in_dim
        self.out_dim = in_dim

        if isinstance(dropout_rate, list):
            assert len(dropout_rate) == len(latent_dim)
        else:
            dropout_rate = [dropout_rate] * len(latent_dim)

        inputs_encoder = tf.keras.Input(shape=(in_dim,), name='counts')
        inputs_sf = tf.keras.Input(shape=(1,), name='size_factors')
        inputs_encoder_pp = PreprocInput()(inputs_encoder)
        output_encoder = Encoder(
            latent_dim=latent_dim[:n_layers_enc],
            l2_coef=l2_coef,
            l1_coef=l1_coef,
            hidden_dropout=dropout_rate[:n_layers_enc],
            input_dropout=input_dropout,
            batchnorm=batchnorm,
            activation=activation,
            kernel_initializer=init
        )(inputs_encoder_pp)

        output_decoder = Decoder(
            latent_dim=latent_dim[n_layers_enc:],
            l2_coef=l2_coef,
            l1_coef=l1_coef,
            hidden_dropout=dropout_rate[n_layers_enc:],
            batchnorm=batchnorm,
            activation=activation,
            kernel_initializer=init
        )(output_encoder)

        if output_layer == 'nb':
            output_decoder_expfamily = NegBinOutput(original_dim=out_dim)((output_decoder, inputs_sf))
        elif output_layer == 'nb_shared_disp':
            output_decoder_expfamily = NegBinSharedDispOutput(original_dim=out_dim)((output_decoder, inputs_sf))
        elif output_layer == 'nb_const_disp':
            output_decoder_expfamily = NegBinConstDispOutput(original_dim=out_dim)((output_decoder, inputs_sf))
        elif output_layer == 'gaussian':
            output_decoder_expfamily = GaussianOutput(original_dim=out_dim)((output_decoder, inputs_sf))
        elif output_layer == 'gaussian_shared_disp':
            output_decoder_expfamily = GaussianSharedStdOutput(original_dim=out_dim)((output_decoder, inputs_sf))
        elif output_layer == 'gaussian_const_disp':
            output_decoder_expfamily = GaussianConstStdOutput(original_dim=out_dim)((output_decoder, inputs_sf))
        else:
            raise ValueError("tried to access a non-supported output layer %s" % output_layer)
        output_decoder_expfamily_concat = tf.keras.layers.Concatenate(axis=1, name="neg_ll")(output_decoder_expfamily)

        self.encoder_model = tf.keras.Model(
            inputs=[inputs_encoder, inputs_sf],
            outputs=output_encoder,
            name="encoder_model"
        )
        self.training_model = tf.keras.Model(
            inputs=[inputs_encoder, inputs_sf],
            outputs=[output_decoder_expfamily_concat],
            name="autoencoder"
        )


class ModelAeVersioned(ModelKerasAe):
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
