from pytest import fixture


def pytest_addoption(parser):
    parser.addoption(
        "--doi_sfaira_repr",
        action="store"
    )
    parser.addoption(
        "--test_data",
        action="store"
    )


@fixture()
def doi_sfaira_repr(request):
    return request.config.getoption("--doi_sfaira_repr")


@fixture()
def test_data(request):
    return request.config.getoption("--test_data")
