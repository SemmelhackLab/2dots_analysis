# Generate the visual prey locations


def prey_location(prey_speed, p1,p2,osci_n):

    prey_speed = prey_speed/300.0 # divide by 0.3 frames per millisecond
    prey = p1
    prey_locations = []
    prey_locations.append(prey)

    if p1 > p2:
        while prey > p2:
            prey -= prey_speed
            prey_locations.append(prey)

    if p1 < p2:
        while prey < p2:
            prey += prey_speed
            prey_locations.append(prey)

    if osci_n > 0:
        if p1 > p2:
            for i in range(osci_n):
                while prey > p2:
                    prey -= prey_speed
                    prey_locations.append(prey)
                while prey < p1:
                    prey += prey_speed
                    prey_locations.append(prey)

        if p1 < p2:
            for i in range(osci_n):
                while prey < p2:
                    prey += prey_speed
                    prey_locations.append(prey)
                while prey > p1:
                    prey -= prey_speed
                    prey_locations.append(prey)

    return prey_locations


def prey_location_osci(prey_speed, p1,p2,p3, dir,osci_n):

    prey_speed = prey_speed/300.0 # divide by 0.3 frames per millisecond
    prey = p1
    prey_locations = []
    prey_locations.append(prey)

    if dir == 'ccw':
        while prey > p2:
            prey -= prey_speed
            prey_locations.append(prey)
        # oscillation
        n = 0
        while n < osci_n:
            while prey < p3:
                prey += prey_speed
                prey_locations.append(prey)
            n += 1
            if n > osci_n-1:
                continue
            while prey > p2:
                prey -= prey_speed
                prey_locations.append(prey)
            n += 1
    if dir == 'cw':
        while prey < p2:
            prey += prey_speed
            prey_locations.append(prey)

        n = 0
        while n < osci_n:
            while prey > p3:
                prey -= prey_speed
                prey_location.append(prey)
            n += 1
            if n > osci_n-1:
                continue
            while prey < p2:
                prey += prey_speed
                prey_locations.append(prey)
            n += 1

    return prey_locations
