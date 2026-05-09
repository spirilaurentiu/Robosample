import robosample


def test_chamber():
    context = robosample.Context(
        name="GfcD",
        seed=42,
        prmtop="examples/GfcDstrippedMin.prmtop",
        inpcrd="examples/GfcDstrippedMin.rst7",
        write_freq=1,
        testing=True,
    )

    # Add cartesian world (will integrate with OpenMM)
    context.add_cartesian_world().add_sampler(
        timeStep=0.001,
        mdSteps=1,
        boostMDSteps=1,
        acceptRejectMode=robosample.rb.AcceptRejectMode.MetropolisHastings,
    )

    sele = context.getDefaultBonds("standard")
    context.add_torsional_world(sele).add_sampler(
        timeStep=0.005, mdSteps=64, boostMDSteps=64
    )

    context.initialize([300.0])
    context.RunREX(0, 1)
