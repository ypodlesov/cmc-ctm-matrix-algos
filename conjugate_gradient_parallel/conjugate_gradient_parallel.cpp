#include "conjugate_gradient_parallel.h"
#include "helpers.h"
#include <thread>
#include <utility>

bool ConjugateGradientParallel(const SparseMatrix<double>& a, const Vector<double>& b, Vector<double>& x) {
    assert(a.data_ && b.data_ && a.row_cnt_ == b.mem_size_);
    if (!x.data_ || x.mem_size_ != a.col_cnt_) {
        x = Vector<double>(a.col_cnt_);
    }
    Vector<double> current_residual(b);
    Vector<double> current_p(b);
    Vector<double> current_x(a.row_cnt_);
    NHelpers::Nullify(current_x.data_, current_x.mem_size_);
    int64_t n = b.mem_size_;

    const int64_t cores_num = std::min<int64_t>(x.mem_size_, std::thread::hardware_concurrency() / 2);
    assert(a.row_cnt_ % cores_num == 0);
    const int64_t row_step = a.row_cnt_ / cores_num;

    double current_alpha, current_beta;

    for (int64_t j = 0; j < x.mem_size_ * x.mem_size_ && !NHelpers::RoughEq<double, double>(current_residual.Norm2(), 0.0); ++j) {
        Vector<double> ap(n);
        {
            std::vector<std::thread> threads;
            threads.reserve(cores_num);
            for (int64_t cur_row_start = 0; cur_row_start < a.row_cnt_; cur_row_start += row_step) {
                threads.emplace_back(&SparseMatrix<double>::VecMult, &a, 
                    std::cref(current_p),
                    std::ref(ap),
                    cur_row_start, 
                    cur_row_start + row_step);
            }
            for (auto&& thread : threads) {
                thread.join();
            }
        }
    
        
        double ap_cur_p_dot_prod = NHelpers::InnerProd(ap.data_, current_p.data_, n);
        if (NHelpers::RoughEq<double, double>(ap_cur_p_dot_prod, 0.0)) {
            break;
        }

        double current_residual_norm = NHelpers::InnerProd(current_residual.data_, current_residual.data_, n);
        current_alpha = current_residual_norm / ap_cur_p_dot_prod;

        current_x.PlusAX(current_p, current_alpha);

        Vector<double> next_residual(n);
        next_residual.AXPlusBY(current_residual, 1, ap, -current_alpha);

        if (NHelpers::RoughEq<double, double>(current_residual_norm, 0.0)) {
            break;
        }

        current_beta = NHelpers::InnerProd(next_residual.data_, next_residual.data_, n) / current_residual_norm;
        current_p.AXPlusBY(next_residual, 1, current_p, current_beta);
        current_residual = std::move(next_residual);
    }
    x = std::move(current_x);
    return true;
}