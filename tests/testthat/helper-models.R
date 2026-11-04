# CHANGE THIS to your package name
PKG <- "IRTnorm"

# Position in this vector = irt_mod number used in data_prep(), fit_IRTnorm(), etc.
ALL_MODELS <- c("1PLnorm", "2PLnorm", "3PLnorm", "GRMnorm",
                "Testlet1PLnorm", "Testlet2PLnorm",
                "HO1PLnorm-C", "HO2PLnorm-C",
                "HO1PLnorm-P", "HO2PLnorm-P",
                "HO1PLnorm-F", "HO2PLnorm-F")

# Position = age_mod number (1 = poly, 2 = splines)
AGE_MODELS <- c("poly", "splines")

MULTI_MODELS <- ALL_MODELS[5:12]