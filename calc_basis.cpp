  // Calculate the basis functions for the given parameter value and knot span
  void CalculateBasisFunctions(int i, float t, std::vector<float>& basis_functions) const {
    // Initialize the left and right arrays
    std::vector<float> left(degree_ + 1);
    std::vector<float> right(degree_ + 1);
    for (int j = 0; j <= degree_; j++) {
      left[j] = t - knots_[i + 1 - j];
      right[j] = knots_[i + j] - t;
      basis_functions[j] = 0.0f;
    }
    basis_functions[0] = 1.0f;
    for (int j = 1; j <= degree_; j++) {
      float saved = 0.0f;
      for (int r = 0; r < j; r++) {
        float temp = basis_functions[r] / (right[r + 1] + left[j - r]);
        basis_functions[r] = saved + right[r + 1] * temp;
        saved = left[j - r] * temp;
      }
      basis_functions[j] = saved;
    }
  }
