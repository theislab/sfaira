Models
======

Introduction to sfaira model management
---------------------------------------

A sfaira model is a class that inherits from BasicModel which defines a tf.keras.models.Model in self.training_model.
This training_model describes the full forward pass. Additionally, embedding models also have an attribute X, a
tf.keras.models.Model that desribes the partial forward pass into the embedding layer.

Such a model class, e.g. ModelX, is wrapped by an inheriting class ModelXVersioned, that handles properties of the
model architecture.
In particular, ModelXVersioned

    - has access to the cell ontology container (a daughter class of CelltypeVersionsBase) that corresponds to this
        model if applicable
    - has access to a map of a version ID to an architectural hyperparameter setting (Topologies), allowing this class
        to set depth, width, etc of the model directly based on the name of the yielded model.
    - has access to the feature space of the model, including its gene names, which are defined by the model topology
        in Topologies


Contribute models
-----------------

Models can be contributed and used in two ways

    - Full model code in sfaira repo
    - Sfaira compatible model code in external package (to come)
