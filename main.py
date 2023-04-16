from parafoil import PFSim


def main():
    lander = PFSim("data/space_rider.yaml")
    lander.start()
    print(lander.state["aerodynamic_force"])
    print(lander.state["aerodynamic_moment"])


if __name__ == "__main__":
    main()
