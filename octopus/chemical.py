NITROUS_OXIDE = 0


def saturation_h(fluid, T, phase):
    if fluid == NITROUS_OXIDE:
        if not (183 < T < 308):
            raise Exception("Temperature must be within range 183-308K")
        Tc = 309.57
        Tr = T / Tc
        if phase == 'l':
            b1 = -200
            b2 = 116.043
            b3 = -917.225
            b4 = 794.779
            b5 = -589.587
        elif phase == 'g':
            b1 = -200
            b2 = 440.055
            b3 = -459.701
            b4 = 434.081
            b5 = -485.338
        else:
            raise Exception("Only phases \'l\' and \'g\' are currently supported")
        return b1 + b2 * (1 - Tr) ** (1 / 3) + b3 * (1 - Tr) ** (2 / 3) + b4 * (1 - Tr) + b5 * (1 - Tr) ** (4 / 3)
