from functools import wraps
from typing import Callable, List


def dataframe_not_empty_check(
    dataframe_arg_index: int = 0,
) -> Callable[[Callable], Callable]:
    """
    A wrapper to validate that a DataFrame is not empty before calling the function.

    :param dataframe_arg_index: The index of the DataFrame argument in the wrapped function's arguments.
    :type dataframe_arg_index: int
    :return: A decorator that wraps the provided function.
    :rtype: Callable[[Callable], Callable]
    """

    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs):
            dataframe = (
                args[dataframe_arg_index]
                if len(args) > dataframe_arg_index
                else kwargs.get("dataframe")
            )
            if dataframe is None:
                raise ValueError("DataFrame argument not provided.")
            if dataframe.empty:
                raise ValueError("Validation failed: DataFrame is empty.")

            return func(*args, **kwargs)

        return wrapper

    return decorator


def dataframe_columns_exist_check(
    columns: List[str], dataframe_arg_index: int = 0
) -> Callable[[Callable], Callable]:
    """
    A wrapper to validate that specific columns exist in a DataFrame before calling the function.

    :param columns: List of column names to validate in the DataFrame.
    :type columns: List[str]
    :param dataframe_arg_index: The index of the DataFrame argument in the wrapped function's arguments.
    :type dataframe_arg_index: int
    :return: A decorator that wraps the provided function.
    :rtype: Callable[[Callable], Callable]
    """

    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs):
            dataframe = (
                args[dataframe_arg_index]
                if len(args) > dataframe_arg_index
                else kwargs.get("dataframe")
            )
            if dataframe is None:
                raise ValueError("DataFrame argument not provided.")

            missing_columns = [col for col in columns if col not in dataframe.columns]
            if missing_columns:
                raise ValueError(
                    f"Validation failed: The following required columns are missing from the DataFrame: {missing_columns}"
                )

            return func(*args, **kwargs)

        return wrapper

    return decorator
