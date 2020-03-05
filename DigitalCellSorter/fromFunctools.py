
class partial:

    """Return "partially evaluated" version of given function/arguments."""

    def __init__(self, func, *args, **kwargs):

        self.func, self.args, self.kwargs = func, args, kwargs

    def __call__(self, *more_args, **more_kwargs):

        all_kwargs = {**self.kwargs, **more_kwargs}

        return self.func(*self.args, *more_args, **all_kwargs)

