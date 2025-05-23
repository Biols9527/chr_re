import yaml

class DefaultConfig:
    """
    Default configuration values for the chr_re package.
    """
    def __init__(self):
        self.reconstruction_methods = ['parsimony', 'ml']
        self.default_reconstruction_method = 'parsimony' # Example: default method to use if multiple are available
        
        # Parsimony specific defaults
        self.parsimony_method_variant = 'fitch' # e.g., 'fitch', 'sankoff'
        
        # Maximum Likelihood specific defaults
        self.ml_model = 'Mk' # e.g., 'Mk', 'BM', 'OU' for continuous
        self.ml_optimization_algorithm = 'L-BFGS-B'
        
        # Bayesian specific defaults
        self.bayesian_mcmc_chains = 2
        self.bayesian_mcmc_generations = 100000
        self.bayesian_mcmc_burnin = 25000
        self.bayesian_prior_type = 'uniform' # e.g., 'uniform', 'exponential'

        # Data processing defaults
        self.missing_data_handling = 'ignore' # e.g., 'ignore', 'infer', 'remove_sites'
        self.data_format = 'nexus' # e.g., 'nexus', 'phylip', 'csv'

        # Analysis defaults
        self.rate_estimation_method = 'mean' # e.g., 'mean', 'median'
        self.event_detection_threshold = 0.05 # p-value or other significance threshold

        # Output and logging
        self.output_directory = './chr_re_output'
        self.log_level = 'INFO' # e.g., 'DEBUG', 'INFO', 'WARNING', 'ERROR'
        self.save_plots = True
        self.plot_format = 'png' # e.g., 'png', 'svg', 'pdf'

    def __str__(self):
        return str(self.__dict__)

def load_config(config_path: str) -> dict:
    """
    Loads a configuration from a YAML file and merges it with default settings.

    Args:
        config_path: Path to the YAML configuration file.

    Returns:
        A dictionary containing the configuration.
    
    Raises:
        FileNotFoundError: If the config_path does not exist.
        yaml.YAMLError: If there is an error parsing the YAML file.
        Exception: For other potential errors during loading.
    """
    # Start with default configuration
    config_data = DefaultConfig().__dict__

    try:
        with open(config_path, 'r') as f:
            user_config = yaml.safe_load(f)
        
        if user_config: # If the user file is not empty
            # Merge user config into default config
            # User's values will overwrite defaults if keys match
            config_data.update(user_config)
            
    except FileNotFoundError:
        print(f"Error: Configuration file not found at '{config_path}'. Using default configuration.")
        # Optionally, re-raise or handle as per application needs
        # raise
    except yaml.YAMLError as e:
        print(f"Error parsing YAML configuration file at '{config_path}': {e}. Using default configuration or last valid state.")
        # Optionally, re-raise or handle
        # raise
    except Exception as e:
        print(f"An unexpected error occurred while loading configuration from '{config_path}': {e}. Using default configuration or last valid state.")
        # Optionally, re-raise or handle
        # raise
        
    return config_data

if __name__ == '__main__':
    # Example usage:
    
    # 1. Create a dummy user config file for testing
    dummy_user_config = {
        'reconstruction_methods': ['ml'],
        'ml_model': 'OU',
        'output_directory': './custom_output',
        'new_user_param': 'test_value'
    }
    dummy_config_path = 'user_config_example.yaml'
    with open(dummy_config_path, 'w') as f:
        yaml.dump(dummy_user_config, f)

    print("--- Loading Default Config (no user file) ---")
    default_settings = DefaultConfig()
    print(default_settings)
    print(f"ML Model (default): {default_settings.ml_model}")

    print("\n--- Loading User Config ---")
    loaded_settings = load_config(dummy_config_path)
    print(loaded_settings)
    print(f"ML Model (user override): {loaded_settings.get('ml_model')}")
    print(f"New User Param: {loaded_settings.get('new_user_param')}")
    print(f"Parsimony Method (default, as not in user file): {loaded_settings.get('parsimony_method_variant')}")

    print("\n--- Testing File Not Found ---")
    non_existent_config = load_config('non_existent_config.yaml')
    # print(non_existent_config) # This will be the default config

    print("\n--- Testing Invalid YAML ---")
    invalid_yaml_path = 'invalid_config.yaml'
    with open(invalid_yaml_path, 'w') as f:
        f.write("reconstruction_methods: ['ml'\nml_model: OU") # Invalid YAML syntax
    
    invalid_settings = load_config(invalid_yaml_path)
    # print(invalid_settings) # This will also be the default config

    # Clean up dummy files
    import os
    os.remove(dummy_config_path)
    os.remove(invalid_yaml_path)

    print("\n--- Verifying DefaultConfig can be accessed directly ---")
    direct_defaults = DefaultConfig()
    print(f"Direct default ML Model: {direct_defaults.ml_model}")
    direct_defaults.ml_model = "custom_test"
    print(f"Modified direct default ML Model: {direct_defaults.ml_model}")
    
    # Ensure load_config returns a dict, not a DefaultConfig instance, but based on it
    print(f"\nType of loaded_settings: {type(loaded_settings)}")
    assert isinstance(loaded_settings, dict)
    assert isinstance(default_settings, DefaultConfig)

    print("\nExample of accessing a config value:")
    current_config = load_config(dummy_config_path) # Re-load for a clean example
    print(f"Current ML Model to be used: {current_config.get('ml_model', default_settings.ml_model)}")

    # Clean up dummy file again if it was recreated by example
    if os.path.exists(dummy_config_path):
        os.remove(dummy_config_path)
    if os.path.exists(invalid_yaml_path):
        os.remove(invalid_yaml_path)
