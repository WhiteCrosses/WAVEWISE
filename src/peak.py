class Peak:
    def __init__(self, frequency, power, distance):
        self.frequency = frequency
        self.distance = distance
        self.min = self.frequency - self.distance
        self.max = self.frequency + self.distance
        self.distanceToNext = None
        self.tickCounter = 0
        self.isChecked = False
        self.power = power
        self.color = 'red'

    def found(self, x, power):
        self.frequency = x
        self.power = power
        self.min = self.frequency - self.distance
        self.max = self.frequency + self.distance
        self.resetTickCounter()

    def resetTickCounter(self):
        self.isChecked = True
        self.tickCounter = 0
