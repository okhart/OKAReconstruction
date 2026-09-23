#ifndef SETTINGS_H
#define SETTINGS_H

// =====================================================
// ===== INPUT / OUTPUT ================================
// =====================================================

/// Name of the input ROOT file containing the event tree
const char* INPUT_FILE = "output2.root";

/// Name of the TTree inside the input file
const char* TREE_NAME  = "decay_points";

/// Number of events to process (0 means all)
const int MAX_EVENTS = 500000;

// =====================================================
// ===== DATA MODE =====================================
// =====================================================

enum class DataMode
{
    Simulation,   ///< Truth information is available for validation
    Experimental  ///< No truth, only reconstruction and output
};

/// Select the data mode
const DataMode DATA_MODE = DataMode::Simulation;

// =====================================================
// ===== TRACK FINDER ==================================
// =====================================================

enum class TrackingAlgorithm
{
    Kalman,  ///< Combinatorial Kalman Filter (recommended)
    Hough    ///< Hough Transform (alternative)
};

/// Choose which tracking algorithm to use
const TrackingAlgorithm TRACKING_ALGORITHM =
    TrackingAlgorithm::Kalman;

// =====================================================
// ===== DEBUG =========================================
// =====================================================

/// Enable verbose debug output and histogram saving for Hough
const bool DEBUG = false;

// =====================================================
// ===== COMBINATORIAL KALMAN FILTER (CKF) =============
// =====================================================

/// Minimum number of hits required to form a track
const int CKF_MIN_HITS = 3;

/// Maximum allowed χ² for a hit to be accepted
const double CKF_CHI2_CUT = 32.0;

/// Process noise added to the slope (simulates multiple scattering)
const double CKF_PROCESS_NOISE = 0.0;

/// Ambiguity threshold: if the difference between the best and second-best
/// χ² is smaller than this, the track is rejected as ambiguous
const double CKF_AMBIGUITY_CHI2 = 0.5;

/// Greedy mode: remove hits after each found track
const bool CKF_GREEDY = true;

/// Generous mode: explore all possible branches (overrides greedy)
const bool CKF_GENEROUS = false;

// =====================================================
// ===== HOUGH TRANSFORM PARAMETERS ====================
// =====================================================

/// Minimum slope (a) in the Hough space
const double HOUGH_A_MIN = -0.09;

/// Maximum slope (a) in the Hough space
const double HOUGH_A_MAX =  0.09;

/// Minimum intercept (b) in the Hough space (if both zero, auto-calc)
const double HOUGH_B_MIN = -80.0;

/// Maximum intercept (b) in the Hough space (if both zero, auto-calc)
const double HOUGH_B_MAX =  80.0;

/// Number of bins in the slope dimension
const int HOUGH_A_BINS = 6;

/// Number of bins in the intercept dimension
const int HOUGH_B_BINS = 13;

/// Fractional overlap between bins (to avoid edge effects)
const double HOUGH_OVERLAP = 0.03;

/// Minimum number of hits in a bin to consider it a valid line candidate
const int HOUGH_THRESHOLD = 6;

// =====================================================
// ===== EVENT QUALITY CUTS ============================
// =====================================================

/// Maximum χ² per degree of freedom for a reconstructed track
const double CHI2_PER_NDF_CUT = 1.0;

/// Expected number of tracks per event (used for validation and filtering)
const int EXPECTED_N_TRACKS = 3;

// =====================================================
// ===== OUTPUT ========================================
// =====================================================

/// Print summary statistics at the end
const bool PRINT_STATISTICS = true;

#endif // SETTINGS_H
