import numpy as np
import os
from typing import Union

from sfaira.ui import UserInterface
from sfaira.unit_tests import DIR_TEMP


class HelperUi:
    ui: Union[UserInterface]
    data: np.ndarray

    """
    Contains functions _test* to test individual functions and attributes of the user ui class.

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

    def test_basic(self):
        """
        Test all relevant model methods.


        :return:
        """
        temp_fn = os.path.join(DIR_TEMP, "test_data")
        self.ui = UserInterface(custom_repo=temp_fn, sfaira_repo=False)


def _test_for_fatal():
    """
    TODO need to simulate/add look up table as part of unit tests locally
    """
    ui = HelperUi()
    ui.test_basic()
