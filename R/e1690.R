#' Phase III cutaneous melanoma clinical trial
#'
#' A well-known dataset from a Phase III cutaneous melanoma clinical trial (Ibrahim et al., 2001).
#'     The study was conducted by the Eastern Cooperative Oncology Group (ECOG) and aimed to evaluate
#'     the effectiveness of the Interferon alpha-2b chemotherapy in preventing recurrence after
#'     surgery. After eliminating missing data, the observations include 417 patients from 1991 to
#'     1995, with follow-up until 1998.
#'
#' @format
#'
#' A data frame with 417 rows and 9 columns:
#' \describe{
#'   \item{trt}{Indicates whether the patient was treated with Interferon alpha-2b chemotherapy
#'   (chemotherapy) or not (control).}
#'   \item{time}{Post-surgery survival or censoring time, in years.}
#'   \item{status}{Censoring status.}
#'   \item{age}{Age, in years.}
#'   \item{sex}{Sex of the patient.}
#'   \item{thickness}{Tumor thickness, in mm.}
#'   \item{nodeII, nodeIII, nodeIV}{Binary variables indicating the nodal category.}
#' }
#'
#' @references
#'
#' Ibrahim, J. G., Chen, M., Sinha, D. (2001). \emph{Bayesian Survival Analysis}. Springer.
#'
#' @examples
#' data(e1690)
#' plot(time ~ trt, e1690, xlab = "Treatment", ylab = "Time")
"e1690"
