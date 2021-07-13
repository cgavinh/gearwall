class PeriodicTable:
    def __init__(self):
        self._numbers = None
        self._symbols = None

    @property
    def numbers(self):
        self._numbers = {1}
        return self._numbers