#pragma once

#include "util.hpp"

/**
 * Computes the symmetric normalized Laplacian of matrix a
 * @param a input matrix
 * @return the symmetric normalized Laplacian
 */
Matd laplacian(const Matd &a);
