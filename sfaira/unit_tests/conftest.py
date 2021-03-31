from pytest import fixture


def pytest_addoption(parser):
    parser.addoption(
        "--doi_sfaira_repr",
        action="store"
    )


@fixture()
def doi_sfaira_repr(request):
    return request.config.getoption("--doi_sfaira_repr")
