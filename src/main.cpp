/*
 * Main executable:
 * Compute the cross sections and print the results to a separate output.
 *
 */

#include "differential_Xsections.cpp"
#include "total_Xsections.cpp"

int main() {
  std::cout << "Checking total xsection:" << std::endl;
  double sigma_C11 = total_xsection_C11(0.994 * 0.994);
  double sigma_C12 = total_xsection_C12(0.994 * 0.994);
  double sigma_C13 = total_xsection_C13(0.994 * 0.994);
  double sigma_C14 = total_xsection_C14(0.994 * 0.994);
  double sigma_C15 = total_xsection_C15(0.994 * 0.994);
  double sigma_C16 = total_xsection_C16(0.994 * 0.994);
  double sigma_C21 = total_xsection_C21(0.994 * 0.994);
  double sigma_C22 = total_xsection_C22(0.994 * 0.994);
  std::cout << sigma_C11 << '\n';
  std::cout << sigma_C12 << '\n';
  std::cout << sigma_C13 << '\n';
  std::cout << sigma_C14 << '\n';
  std::cout << sigma_C15 << '\n';
  std::cout << sigma_C16 << '\n';
  std::cout << sigma_C21 << '\n';
  std::cout << sigma_C22 << '\n';
  std::cout << "Check :)" << std::endl;
  std::cout << " " << std::endl;

  std::cout << "Checking differential xsection:" << std::endl;
  double diff_sigma_C11 = diff_xsection_C11(-0.25, 0.994 * 0.994);
  double diff_sigma_C12 = diff_xsection_C12(-0.25, 0.994 * 0.994);
  double diff_sigma_C13 = diff_xsection_C13(-0.25, 0.994 * 0.994);
  double diff_sigma_C14 = diff_xsection_C14(-0.25, 0.994 * 0.994);
  double diff_sigma_C15 = diff_xsection_C15(-0.25, 0.994 * 0.994);
  double diff_sigma_C16 = diff_xsection_C16(-0.25, 0.994 * 0.994);
  double diff_sigma_C21 = diff_xsection_C21(-0.25, 0.994 * 0.994);
  double diff_sigma_C22 = diff_xsection_C22(-0.25, 0.994 * 0.994);
  std::cout << diff_sigma_C11 << '\n';
  std::cout << diff_sigma_C12 << '\n';
  std::cout << diff_sigma_C13 << '\n';
  std::cout << diff_sigma_C14 << '\n';
  std::cout << diff_sigma_C15 << '\n';
  std::cout << diff_sigma_C16 << '\n';
  std::cout << diff_sigma_C21 << '\n';
  std::cout << diff_sigma_C22 << '\n';
  std::cout << "Check :)" << std::endl;
}
