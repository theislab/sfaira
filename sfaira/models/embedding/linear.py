import numpy as np
import tensorflow as tf
from typing import List, Union

from sfaira.models.embedding.output_layers import NegBinOutput, NegBinSharedDispOutput, NegBinConstDispOutput, \
    GaussianOutput, GaussianSharedStdOutput, GaussianConstStdOutput
from sfaira.models.embedding.external import BasicModel
from sfaira.models.embedding.external import PreprocInput
from sfaira.models.embedding.external import Topologies


class EncoderLinear(tf.keras.layers.Layer):
    """Maps input to embedding space"""

    def __init__(
            self,
            latent_dim: int,
            positive_components: bool,
            l1_coef: float,
            l2_coef: float,
            name='encoder',
            **kwargs
    ):
        super().__init__(name=name, **kwargs)
        self.fwd_pass = []
        self.fwd_pass.append(tf.keras.layers.Dense(
            units=latent_dim,
            activation="linear" if not positive_components else "relu",
            kernel_regularizer=tf.keras.regularizers.l1_l2(l1=l1_coef, l2=l2_coef)
        ))

    def call(self, inputs, **kwargs):
        x = inputs
        for layer in self.fwd_pass:
            x = layer(x)
        return x


class ModelLinear(BasicModel):

    def __init__(
            self,
            in_dim,
            latent_dim: int = 10,
            positive_components: bool = False,
            l2_coef: float = 0.,
            l1_coef: float = 0.,
            dropout_rate=None,
            output_layer="nb"
    ):
        super(ModelLinear, self).__init__()

        self.in_dim = in_dim
        self.latent_dim = latent_dim
        out_dim = in_dim
        self.out_dim = in_dim

        inputs_encoder = tf.keras.Input(shape=(in_dim,), name='counts')
        inputs_sf = tf.keras.Input(shape=(1,), name='size_factors')
        inputs_encoder_pp = PreprocInput()(inputs_encoder)
        output_encoder = EncoderLinear(
            latent_dim=latent_dim,
            positive_components=positive_components,
            l1_coef=l1_coef,
            l2_coef=l2_coef
        )(inputs_encoder_pp)

        if output_layer == 'nb':
            output_decoder_expfamily = NegBinOutput(original_dim=out_dim)((output_encoder, inputs_sf))
        elif output_layer == 'nb_shared_disp':
            output_decoder_expfamily = NegBinSharedDispOutput(original_dim=out_dim)((output_encoder, inputs_sf))
        elif output_layer == 'nb_const_disp':
            output_decoder_expfamily = NegBinConstDispOutput(original_dim=out_dim)((output_encoder, inputs_sf))
        elif output_layer == 'gaussian':
            output_decoder_expfamily = GaussianOutput(original_dim=out_dim)((output_encoder, inputs_sf))
        elif output_layer == 'gaussian_shared_disp':
            output_decoder_expfamily = GaussianSharedStdOutput(original_dim=out_dim)((output_encoder, inputs_sf))
        elif output_layer == 'gaussian_const_disp':
            output_decoder_expfamily = GaussianConstStdOutput(original_dim=out_dim)((output_encoder, inputs_sf))
        else:
            raise ValueError("tried to access a non-supported output layer %s" % output_layer)
        output_decoder_expfamily_concat = tf.keras.layers.Concatenate(axis=1, name="neg_ll")(output_decoder_expfamily)

        self.encoder_model = tf.keras.Model(
            inputs=inputs_encoder,
            outputs=output_encoder,
            name="encoder"
        )
        self.training_model = tf.keras.Model(
            inputs=[inputs_encoder, inputs_sf],
            outputs=output_decoder_expfamily_concat,
            name="autoencoder"
        )

    def predict_reconstructed(self, x: np.ndarray):
        return np.split(self.training_model.predict(x), indices_or_sections=2, axis=1)[0]

    def predict_embedding(self, x: np.ndarray, **kwargs):
        return self.encoder.predict(x)


class ModelLinearVersioned(ModelLinear):
    def __init__(
            self,
            topology_container: Topologies,
            override_hyperpar: Union[dict, None] = None
    ):
        hyperpar = topology_container.topology["hyper_parameters"]
        if override_hyperpar is not None:
            for k in list(override_hyperpar.keys()):
                hyperpar[k] = override_hyperpar[k]
        ModelLinear.__init__(
            self=self,
            in_dim=topology_container.ngenes,
            **hyperpar
        )
        print('passed hyperpar: \n', hyperpar)
        self._topology_id = topology_container.topology_id
        self.genome_size = topology_container.ngenes
        self.model_class = topology_container.model_class
        self.model_type = topology_container.model_type
        self.hyperparam = dict(
            list(hyperpar.items()) +
            [
                ("topology_id", self._topology_id),
                ("genome_size", self.genome_size),
                ("model_class", self.model_class),
                ("model_type", self.model_type)
            ]
        )
