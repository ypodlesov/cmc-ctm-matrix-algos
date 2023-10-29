class TMultiplyTestBase : public testing::TestWithParam<size_t> {
public:
    TMultiplyTestBase() 
        : N{GetParam()}
    {
    }

    virtual void SetUp() override {
        std::srand(std::time(nullptr));
        A = TMatrix<double>::CreateRandom(N, N);
        B = TMatrix<double>::CreateRandom(N, N);
    } 

protected:
    const size_t N;
    TMatrix<double> A;
    TMatrix<double> B;
};