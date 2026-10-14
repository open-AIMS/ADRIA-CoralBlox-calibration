using ADRIACoralBloxCalibration
using ADRIACoralBloxCalibration.common
using ADRIACoralBloxCalibration.calibration

config = load_config(joinpath(@__DIR__, "..", "config.toml"))
dom = load_domain(config)
location_classification = load_location_classification(config.loc_class_path)

cfg = CalibConfig(dom)

calib_data = build_calibration_data(
    dom, config.ltmp_reef_data_path, config.composition_path;
    out_dir=config.out_dir
)

run_calibration(
    dom,
    cfg,
    calib_data,
    location_classification.consecutive_classification,
    10_000;
    init_cover_path=config.init_cover_path,
    out_dir=config.out_dir,
    init_guess_path=config.init_guess_path,
    config=config,
)
