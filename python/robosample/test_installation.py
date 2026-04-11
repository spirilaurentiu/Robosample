from pathlib import Path


def run_test():
    import robosample

    current_dir = Path(__file__).parent.resolve()
    repo_root = current_dir.parents[1]
    examples_dir = repo_root / "examples"

    # Create simulation
    name = "1APQ.test.rigid"
    seed = 42
    prmtop = str(examples_dir / "1APQ.prmtop")
    xyz = str(examples_dir / "1APQ.rst7")
    write_freq = 0
    equil_steps = 0
    prod_steps = 0
    context = robosample.Context(
        name=name,
        seed=seed,
        prmtop=prmtop,
        inpcrd=xyz,
        write_freq=write_freq,
        testing=True,
    )

    # Add cartesian world (will integrate with OpenMM)
    context.addCartesianWorld().addSampler(
        timeStep=0.001,
        mdSteps=10,
        boostMDSteps=10,
        acceptRejectMode=robosample.rb.AcceptRejectMode.AlwaysAccept,
    )

    # Add torsional world
    sele = context.getDefaultBonds("standard")
    context.addTorsionalWorld(sele).addSampler(
        timeStep=0.01, mdSteps=10, boostMDSteps=10
    )

    # Add one replica at 300 K
    context.initialize([300.0])

    # Run the simulation for no equilibration steps and 100 production steps totaling 100 ps of simulation time
    context.RunREX(equil_steps, prod_steps)


if __name__ == "__main__":
    try:
        run_test()
    except Exception as e:
        print(f"Error: {e}")
    else:
        print("\n\nRobosample installed successfully.")
