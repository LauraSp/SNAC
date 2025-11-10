import argparse

from snac.SNACmodel import AggregationModel
from snac.diamond import Diamond


def model(save_dir: str,
          model_data=None,
          diamond_data=None, cooling_rate0=None, T_start0=None,
          rate_bounds=None, T_bounds=None,
          ):
    """Run AggregationModel with specified parameters.

    Either a diamond_file and cooling history parameters are provided, or a
    previously saved model file is loaded.

    PARAMS:
    --------
    """
    # load model or create a new one depending on what the input is
    if model_data is not None:
        if isinstance(model_data, str):
            # Load model from file
            model = AggregationModel.from_json(model_data)
        elif isinstance(model_data, AggregationModel):
            model = model_data
        else:
            raise ValueError(
                "model_data must be a file path or AggregationModel instance."
                )

    else:
        # Load diamond data from file
        if isinstance(diamond_data, str):
            diamond = Diamond.from_json(diamond_data)

        elif isinstance(diamond_data, Diamond):
            diamond = diamond_data
        else:
            raise ValueError(
                "diamond_data must be a file path or Diamond instance."
                )

        if not all(v is not None for v in
                   [cooling_rate0, T_start0, rate_bounds, T_bounds]):
            raise ValueError(
                "cooling_rate0, T_start0, rate_bounds, and T_bounds "
                "must be provided when diamond_data is given."
                )

        # Create AggregationModel instance
        model = AggregationModel(
            diamond=diamond,
            cooling_rate0=cooling_rate0,
            T_start0=T_start0,
            rate_bounds=rate_bounds,
            T_bounds=T_bounds
         )

    # Run the model
    model.optimise_history()

    # Print results
    print(
        f"Model result:\ninitial T: {model.model_results[0]:.0f} deg.C\n"
        f"cooling rate: {1000*model.model_results[1]:.0f} K/Gyr")

    # Save results
    model.to_json(save_dir + "/model_results.json")
    model.save_history(save_dir + "/model_history.csv")

    print(f"Model results saved to {save_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
       "--save_dir",
       type=str,
       required=True,
       help="Directory to save model results"
       )
    parser.add_argument(
       "--model_file",
       type=str,
       required=False,
       help=("Path to the model file. If provided, "
             "other parameters are ignored.")
    )
    parser.add_argument(
       "--diamond_file",
       type=str,
       required=False,
       help="Path to the diamond file"
       )
    parser.add_argument(
       "--cooling_rate0",
       type=float,
       required=False,
       help="Initial cooling rate (K/Gyr)"
       )
    parser.add_argument(
       "--T_start0",
       type=float,
       required=False,
       help="Initial temperature (deg.C)"
       )
    parser.add_argument(
       "--rate_bounds",
       type=float,
       nargs=2,
       required=False,
       help="Cooling rate bounds (min, max)"
       )
    parser.add_argument(
       "--T_bounds",
       type=float,
       nargs=2,
       required=False,
       help="Temperature bounds (min, max)"
       )

    args = parser.parse_args()
    print(args)

    model(
        save_dir=args.save_dir,
        model_data=args.model_file,
        diamond_data=args.diamond_file,
        cooling_rate0=args.cooling_rate0,
        T_start0=args.T_start0,
        rate_bounds=args.rate_bounds,
        T_bounds=args.T_bounds
    )
