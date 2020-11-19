from choose_sequences import sequences

sequences("xenonbaseline_highres.csv", False, 1, 1, 150, False, True, True)
sequences("argon.csv", "xenonbaseline_highres.csv", 1, 0, 150, False, True, True)
sequences("mix.csv", "xenonbaseline_highres.csv", 1, 0, 150, False, True, True)
sequences("xenon.csv", "xenonbaseline_highres.csv", 1, 0, 150, False, True, True)
