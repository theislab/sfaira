
def crossref_query(doi: str, k: str):
    """
    Queries cross REST API via package crossref_commons.

    :param doi: DOI to query from crossref for.
    :param k: Key to extract from crossref query container.
    :return:
    """
    from crossref_commons.retrieval import get_entity
    from crossref_commons.types import EntityType, OutputType

    try:
        attempt_counter = 0
        x = None
        while x is None:
            try:
                attempt_counter += 1
                x = get_entity(doi, EntityType.PUBLICATION, OutputType.JSON)[k]
            except ConnectionError as e:
                # Terminate trial after 5 attempts with ConnectionError:
                if attempt_counter > 5:
                    raise ConnectionError(e)
        return x
    except ValueError as e:
        print(f"ValueError: {e}")
        return None
    except ConnectionError as e:
        print(f"ConnectionError: {e}")
        return None
