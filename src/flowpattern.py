from enum import Enum


class FlowPattern(Enum):
    distributed = 1
    intermittent = 2
    transition = 3
    segregated = 4
    downward = 5
