#ifdef USE_GSL

#include "regression.h"
#include <gsl/gsl_multifit.h>

using namespace std;

vector<double> compute_linear_regression(const vector<WeightedValue> &m_data, const vector<const vector<double>*> &m_features,
    const bool m_use_bias)
{
    if( m_features.empty() ){
        throw __FILE__ ":compute_linear_regression: No feature vectors found";
    }

    const size_t num_data = m_data.size();

    if(num_data != m_features.size()){
        throw __FILE__ ":compute_linear_regression: |data| != |feature vectors|";
    }

    size_t num_features = 0;

    for(vector< const vector<double>* >::const_iterator i = m_features.begin();i != m_features.end();++i){

        if((*i) == NULL){
            throw __FILE__ ":compute_linear_regression: Yiikes! A null pointer in the feature vectors ...";
        }

        if( i == m_features.begin() ){

            num_features = (*i)->size();

            if(num_features == 0){
                throw __FILE__ ":compute_linear_regression: Empty feature vector";
            }
        }

        if( (*i)->size() != num_features ){
            throw __FILE__ ":compute_linear_regression: Variable number of features is not allowed";
        }
    }

    // Use the Gnu Scientific Library to compute the linear regression
    // See: https://www.gnu.org/software/gsl/doc/html/lls.html#linear-regression

    const size_t total_num_features = num_features + (m_use_bias ? 1 /* include the bias */ : 0);

    double chisq;
    
    gsl_matrix *X = gsl_matrix_alloc(num_data, total_num_features);
    gsl_vector *y = gsl_vector_alloc(num_data);
    gsl_vector *w = gsl_vector_alloc(num_data);

    gsl_vector *c = gsl_vector_alloc(total_num_features);
    gsl_matrix *cov = gsl_matrix_alloc(total_num_features, total_num_features);

    const size_t bias_offset = (m_use_bias ? 1 : 0);

    // Pack the data into the GSL data structures
    for(size_t i = 0;i < num_data;++i){

        if(m_use_bias){
            gsl_matrix_set(X, i, 0, 1.0); // The bias
        }

        for(size_t j = 0;j < num_features;++j){
            gsl_matrix_set(X, i, j + bias_offset, (* (m_features[i]) )[j]);
        }

        gsl_vector_set(y, i, m_data[i].value);
        gsl_vector_set(w, i, m_data[i].weight);
    }

    {
        gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc(num_data, total_num_features);
        
        gsl_multifit_wlinear(X, w, y, c, cov, &chisq, work);

        gsl_multifit_linear_free (work);
    }

    vector<double> ret(total_num_features);

    // Convert from a GSL vector to a C++ vector
    for(size_t i = 0;i < total_num_features;++i){
        ret[i] = gsl_vector_get(c ,i);
    }

    gsl_matrix_free(X);
    gsl_vector_free(y);
    gsl_vector_free(w);
    gsl_vector_free(c);
    gsl_matrix_free(cov);
    
    return ret;
}

double linear_prediction(const vector<double> &m_model, const vector<double> &m_features)
{
    const bool has_bias = ( m_model.size() == (m_features.size() + 1) );

    if( !has_bias && (m_model.size() != m_features.size() ) ){
        throw __FILE__ ":linear_prediction: |linear parameters| != |features|";
    }

    double ret = has_bias ? m_model[0] : 0.0;
    const size_t num_features = m_features.size();
    size_t bias_offset = has_bias ? 1 : 0;

    for(size_t i = 0;i < num_features;++i){
        ret += m_model[i + bias_offset]*m_features[i];
    }

    return ret;
}

#endif // USE_GSL