# Generates data/calibration_split.csv from a bare ADRIA domain, not via
# ADRIACoralBloxCalibration.load_domain(config) - that attaches historical_disturbance_mortality_rates.nc,
# the very file this pipeline builds, which would be circular on a fresh checkout.
# build_calibration_data only reads dom.loc_data, so a bare domain is sufficient.
#
# Included (guarded on data/calibration_split.csv's existence) from base.jl - relies on
# RME_DOMAIN_PATH, LTMP_REEF_DATA_PATH, COMPOSITION_PATH and DATA_DIR being defined there.
# Delete data/calibration_split.csv to force regeneration.

using ADRIACoralBloxCalibration
using ADRIACoralBloxCalibration.common: START_YEAR, END_YEAR
using ADRIACoralBloxCalibration.calibration: build_calibration_data

@info "No local calibration_split.csv found - generating one from a bare RME domain"
bare_domain = ADRIA.load_domain(
    ADRIA.RMEDomain, RME_DOMAIN_PATH, "45", timeframe=(START_YEAR, END_YEAR)
)
build_calibration_data(
    bare_domain, LTMP_REEF_DATA_PATH, COMPOSITION_PATH; out_dir=DATA_DIR
)
