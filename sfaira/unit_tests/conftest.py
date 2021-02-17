"""
Parameterizing test according to https://stackoverflow.com/questions/40880259/how-to-pass-arguments-in-pytest-by-command-line
"""


def pytest_addoption(parser):
    parser.addoption("--doi_sfaira_repr", action="store", default="d10_1016_j_cmet_2019_01_021")


def pytest_generate_tests(metafunc):
    # This is called for every test. Only get/set command line arguments
    # if the argument is specified in the list of test "fixturenames".
    option_value = metafunc.config.option.name
    if "doi_sfaira_repr" in metafunc.fixturenames and option_value is not None:
        metafunc.parametrize("doi_sfaira_repr", [option_value])
