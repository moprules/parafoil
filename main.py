from parafoil import PFSim


def main():
    lander = PFSim("data/space_rider.yaml")
    lander.start()
    print(lander.state["force_payload_bodyframe"])
    print(lander.state["moment_payload_bodyframe"])


if __name__ == "__main__":
    main()
