void export_2D() {
    TFile *f = TFile::Open("K47_T1_002.root");
    TH2 *h2 = (TH2*) f->Get("hh_gE_ch_beta");

    if (!h2) {
        std::cout << "Histogram not found!" << std::endl;
        return;
    }

    std::ofstream out("hist2d_export.txt");

    int nx = h2->GetNbinsX();
    int ny = h2->GetNbinsY();

    out << "# x_center y_center content\n";

    for (int ix = 1; ix <= nx; ix++) {
        double x = h2->GetXaxis()->GetBinCenter(ix);
        for (int iy = 1; iy <= ny; iy++) {
            double y = h2->GetYaxis()->GetBinCenter(iy);

            double content = h2->GetBinContent(ix, iy);

            out << x << "  " << y << "  " << content << "\n";
        }
    }

    out.close();
    std::cout << "Saved to hist2d_export.txt\n";
}