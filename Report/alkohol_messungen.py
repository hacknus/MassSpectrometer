from choose_sequences import sequences

sequences("ethanol.csv",'vodka.csv', 5, 0, 150, False, True)
sequences("vodka.csv",'air_background.csv', 5, 0, 150, False, True)
sequences("vodka_pure.csv",'air_background.csv', 5, 0, 150, False, True)
#sequences("Xenon.csv",'air_background', 5, 0, 150, False, True)
