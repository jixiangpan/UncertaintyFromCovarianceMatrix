void test_matrix()
{
  TString roostr = "";
  
  gStyle->SetOptStat(0);

  double array_sigma[2] = {2,5};
  
  TMatrix matrix_cor(2,2);
  matrix_cor(0,0) = 1;
  matrix_cor(1,1) = 1;
  matrix_cor(0,1) = 0.2;
  matrix_cor(1,0) = 0.2;
  
  TMatrixD matrix_cov(2,2);
  TMatrixDSym DSmatrix_cov(2);    
  for(int i=0; i<2; i++) {
    for(int j=0; j<2; j++) {
      matrix_cov(i,j) = matrix_cor(i,j) * array_sigma[i] * array_sigma[j];
      DSmatrix_cov(i,j) = matrix_cov(i,j);
    }
  }

  TMatrixDSymEigen DSmatrix_eigen(DSmatrix_cov);
  TMatrixD matrix_eigenvector = DSmatrix_eigen.GetEigenVectors();
  TMatrixD matrix_eigenvector_T(matrix_eigenvector.GetNrows(), matrix_eigenvector.GetNcols());
  matrix_eigenvector_T.Transpose( matrix_eigenvector);
  TMatrixD matrix_diag = matrix_eigenvector_T * matrix_cov * matrix_eigenvector;
  for(int i=0; i<2; i++) {
    for(int j=0; j<2; j++) {
      if( matrix_diag(i,j)<1e-8 ) matrix_diag(i,j) = 0;
    }
  }
  TMatrixD matrix_cov_AA = matrix_eigenvector * matrix_diag * matrix_eigenvector_T;

  TPrincipal p_test(2, "ND");
  double array_data[2] = {0};

  TRandom3 *rand = new TRandom3();
  
  for(int i=1; i<=500000; i++) {
    double x1 = rand->Gaus( 0, sqrt( matrix_diag(0,0) ) );
    double x2 = rand->Gaus( 0, sqrt( matrix_diag(1,1) ) );

    TMatrixD matrix_ele(2,1);
    matrix_ele(0,0) = x1;
    matrix_ele(1,0) = x2;

    TMatrixD matrix_flu = matrix_eigenvector * matrix_ele;

    array_data[0] = matrix_flu(0,0);
    array_data[1] = matrix_flu(1,0);

    // array_data[0] = x1;
    // array_data[1] = x2;

    p_test.AddRow( array_data );
  }
  
  TMatrixD *matrix_cov_BB = p_test.GetCovarianceMatrix();
  
  
  //////
  roostr = "matrix_cor";
  TCanvas *canv_matrix_cor = new TCanvas(roostr, roostr, 900*8./9, 650*8./9);
  canv_matrix_cor->SetRightMargin(0.15);
  matrix_cor.Draw("colz text");

  roostr = "matrix_cov";
  TCanvas *canv_matrix_cov = new TCanvas(roostr, roostr, 900*8./9, 650*8./9);
  canv_matrix_cov->SetRightMargin(0.15);
  matrix_cov.Draw("colz text");

  roostr = "matrix_diag";
  TCanvas *canv_matrix_diag = new TCanvas(roostr, roostr, 900*8./9, 650*8./9);
  canv_matrix_diag->SetRightMargin(0.15);
  matrix_diag.Draw("colz text");

  roostr = "matrix_cov_AA";
  TCanvas *canv_matrix_cov_AA = new TCanvas(roostr, roostr, 900*8./9, 650*8./9);
  canv_matrix_cov_AA->SetRightMargin(0.15);
  matrix_cov_AA.Draw("colz text");

  roostr = "matrix_cov_BB";
  TCanvas *canv_matrix_cov_BB = new TCanvas(roostr, roostr, 900*8./9, 650*8./9);
  canv_matrix_cov_BB->SetRightMargin(0.15);
  matrix_cov_BB->Draw("colz text");



  
}
