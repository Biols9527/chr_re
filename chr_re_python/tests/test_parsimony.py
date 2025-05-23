import pytest
import numpy as np
from chr_re.methods.parsimony import ParsimonyReconstruction
from ete3 import Tree # Assuming ete3 is used for tree loading
import pandas as pd # Assuming pandas for data loading
import os

# Define base path for examples, assuming tests are run from project root or tests/
# This gets to the parent directory of 'tests/', which should be 'chr_re_python'
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLES_DIR = os.path.join(BASE_DIR, "examples")

class TestParsimonyMethods:
    """
    Tests for parsimony-based ancestral state reconstruction methods.
    """
    def setup_method(self, method):
        """
        Set up test data and reconstructor instance before each test method.
        """
        tree_path = os.path.join(EXAMPLES_DIR, "simulated_tree.nwk")
        counts_path = os.path.join(EXAMPLES_DIR, "simulated_counts.csv")

        if not os.path.exists(tree_path):
            raise FileNotFoundError(f"Test setup failed: Tree file not found at {tree_path}")
        if not os.path.exists(counts_path):
            raise FileNotFoundError(f"Test setup failed: Counts file not found at {counts_path}")

        # Load tree (using ete3)
        self.tree = Tree(tree_path)
        # Load counts (using pandas)
        counts_df = pd.read_csv(counts_path)
        self.tip_states = pd.Series(counts_df.chromosome_number.values, index=counts_df.species).to_dict()
        
        # Initialize reconstructor
        self.reconstructor = ParsimonyReconstruction()

    def test_fitch_algorithm_placeholder(self):
        """
        Test that Fitch algorithm raises NotImplementedError.
        """
        with pytest.raises(NotImplementedError):
            self.reconstructor.fitch_algorithm(self.tree, self.tip_states)

    def test_wagner_parsimony_placeholder(self):
        """
        Test that Wagner parsimony raises NotImplementedError.
        """
        with pytest.raises(NotImplementedError):
            self.reconstructor.wagner_parsimony(self.tree, self.tip_states)

    def test_sankoff_algorithm_placeholder(self):
        """
        Test that Sankoff algorithm raises NotImplementedError or ValueError.
        """
        # Scenario 1: No cost matrix provided at all (should raise ValueError in actual method)
        # For placeholder, we test the NotImplementedError when a cost matrix is supplied,
        # as per the example's focus. The ValueError for missing matrix is inherent to Sankoff's logic.
        # If the method only raises NotImplementedError, we test that.
        # The ParsimonyReconstruction.sankoff_algorithm raises ValueError if cost_matrix is None,
        # otherwise it raises NotImplementedError.

        # Test ValueError if no cost matrix is provided
        with pytest.raises(ValueError): # Exact match message can be added if needed
            self.reconstructor.sankoff_algorithm(self.tree, self.tip_states)

        # Test NotImplementedError if a cost matrix is provided (as per placeholder logic)
        dummy_cost_matrix = np.array([[0, 1], [1, 0]]) # Example
        reconstructor_with_matrix_init = ParsimonyReconstruction(cost_matrix=dummy_cost_matrix)
        with pytest.raises(NotImplementedError):
            reconstructor_with_matrix_init.sankoff_algorithm(self.tree, self.tip_states)
        
        # Also test when cost_matrix is passed directly to the method
        with pytest.raises(NotImplementedError):
            self.reconstructor.sankoff_algorithm(self.tree, self.tip_states, cost_matrix=dummy_cost_matrix)


# Example of a fixture if complex setup is needed for multiple tests
# @pytest.fixture
# def sample_tree_and_states():
#     tree_path = os.path.join(EXAMPLES_DIR, "simulated_tree.nwk")
#     counts_path = os.path.join(EXAMPLES_DIR, "simulated_counts.csv")
#     tree = Tree(tree_path)
#     counts_df = pd.read_csv(counts_path)
#     tip_states = pd.Series(counts_df.chromosome_number.values, index=counts_df.species).to_dict()
#     return tree, tip_states

# Example of a test using the fixture
# def test_fitch_with_fixture(sample_tree_and_states):
#     tree, tip_states = sample_tree_and_states
#     reconstructor = ParsimonyReconstruction()
#     with pytest.raises(NotImplementedError):
#         reconstructor.fitch_algorithm(tree, tip_states)
