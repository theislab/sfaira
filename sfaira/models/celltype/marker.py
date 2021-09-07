try:
    import tensorflow as tf
except ImportError:
    tf = None
from typing import Union

from sfaira.versions.metadata import CelltypeUniverse
from sfaira.versions.topologies import TopologyContainer
from sfaira.models.celltype.base import BasicModelKerasCelltype
from sfaira.models.pp_layer import PreprocInput


class LearnedThresholdLayer(tf.keras.layers.Layer):
    """
    A layer that thresholds the input with a learned threshold.
    """

    def __init__(
        self,
        out_dim,
        initializer="zeros",
        name='learned_threshold_layer',
        **kwargs
    ):
        super().__init__(name=name, **kwargs)
        self.bias = self.add_weight(
            "bias",
            shape=[1, out_dim],
            initializer=initializer
        )
        self.slope = self.add_weight(
            "slope",
            shape=[1, out_dim],
            initializer="zeros"
        )

    def call(self, inputs):
        x = inputs
        x = (x + self.bias) * tf.exp(self.slope)
        return tf.nn.sigmoid(x)


class CellTypeMarker(BasicModelKerasCelltype):
    """
    Marker gene-based cell type classifier: Learns whether or not each gene exceeds requires threshold
    and learns cell type assignment as linear combination of these marker gene presence probabilities.
    Activitiy and weight regularizers keep this sparse.
    """

    def __init__(
            self,
            in_dim: int,
            out_dim: int,
            l1_coef: float = 0.,
            l2_coef: float = 0.,
            kernel_initializer='glorot_uniform',
            bias_initializer='zeros',
            bias_regularizer=None,
            kernel_constraint=None,
            bias_constraint=None,
            name='celltype_marker',
            **kwargs
    ):
        """

        :param in_dim: Input feature dimension (number of features).
        :param out_dim: Number of output nodes: Number of cell types +1 for unkown cell type.
        :param l1_coef: L1 regularization hyperparameter on weights of last layer.
        :param l2_coef: L2 regularization hyperparameter on weights of last layer and activation of hidden.
        :param kernel_initializer: See tf.keras.layers.Dense() documentation.
        :param bias_initializer: See tf.keras.layers.Dense() documentation.
        :param kernel_regularizer: See tf.keras.layers.Dense() documentation.
        :param bias_regularizer: See tf.keras.layers.Dense() documentation.
        :param kernel_constraint: See tf.keras.layers.Dense() documentation.
        :param bias_constraint: See tf.keras.layers.Dense() documentation.
        """
        self.in_dim = in_dim
        self.out_dim = out_dim

        inputs = tf.keras.Input(shape=(in_dim,))
        x = PreprocInput()(inputs)
        # Expand one dimension to be able to use LocallyConnected1D.
        x = LearnedThresholdLayer(out_dim=self.in_dim)(x)
        y = tf.keras.layers.Dense(
            units=out_dim,
            activation="softmax",
            use_bias=True,
            kernel_initializer=kernel_initializer,
            bias_initializer=bias_initializer,
            kernel_regularizer=tf.keras.regularizers.l1_l2(l1=l1_coef, l2=l2_coef),
            bias_regularizer=bias_regularizer,
            activity_regularizer=None,
            kernel_constraint=kernel_constraint,
            bias_constraint=bias_constraint
        )(x)
        self.training_model = tf.keras.Model(inputs=inputs, outputs=y, name=name)


class CellTypeMarkerVersioned(CellTypeMarker):

    def __init__(
            self,
            celltypes_version: CelltypeUniverse,
            topology_container: TopologyContainer,
            override_hyperpar: Union[dict, None] = None
    ):
        """

        :param celltypes_version:
        :param topology_container:
        :param override_hyperpar: Dictionary with hyper-parameters of model to override in preset hyper-parameter
            dictionary that is queried based on the topology_id. Can contain a subset of all hyperparameters.
        """
        # Get cell type version instance based on topology ID, organism and organ.
        hyperpar = topology_container.topology["hyper_parameters"]
        if override_hyperpar is not None:
            for k in list(override_hyperpar.keys()):
                hyperpar[k] = override_hyperpar[k]
        super().__init__(
            in_dim=topology_container.n_var,
            out_dim=celltypes_version.onto_cl.n_leaves,
            **hyperpar
        )
        print('passed hyperpar: \n', hyperpar)
        self._topology_id = topology_container.topology_id
        self.genome_size = topology_container.n_var
        self.model_class = "celltype"
        self.model_type = topology_container.model_type
        self.hyperparam = dict(
            list(hyperpar.items()) +  # noqa: W504
            [
                ("topology_id", self._topology_id),
                ("genome_size", self.genome_size),
                ("model_class", self.model_class),
                ("model_type", self.model_type),
                ("ntypes", celltypes_version.onto_cl.n_leaves),
            ]
        )
