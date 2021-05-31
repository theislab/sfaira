import numpy as np
import os
from typing import Union

from sfaira.interface import UserInterface


class TestUi:
    ui: Union[UserInterface]
    data: np.ndarray

    """
    Contains functions _test* to test individual functions and attributes of the user interface class.

    TODO for everybody working on this, add one _test* function in here and add it into
    basic_estimator_test(). See _test_call() for an example.
    """

    def simulate(self):
        """
        Simulate basic data example used for unit test.

        Sets attribute .data with simulated data.

        :return:
        """
        pass

    def _test_basic(self):
        """
        Test all relevant model methods.


        :return:
        """
        temp_fn = os.path.join(str(os.path.dirname(os.path.abspath(__file__))), '../test_data')
        self.ui = UserInterface(custom_repo=temp_fn, sfaira_repo=False)
