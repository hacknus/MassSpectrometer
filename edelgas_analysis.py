from choose_sequences import sequences

sequences("xenonbaseline.csv", False, 5, 1, 150, False, True)
sequences("Argon.csv", "xenonbaseline.csv", 5,0, 150, False, True)
sequences("gemisch.csv", "xenonbaseline.csv", 5, 0, 150, False, True)
sequences("Xenon.csv", "xenonbaseline.csv", 5, 0, 150, False, True)
