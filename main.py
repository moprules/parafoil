from parafoil import PFSim


def main():
    lander = PFSim("data/space_rider.yaml")
    lander.start()
    print(lander.state["mesh"]["Vinf"])


if __name__ == "__main__":
    main()
