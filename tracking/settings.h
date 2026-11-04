#ifndef SETTINGS_H
#define SETTINGS_H

// =====================================================
// ===== INPUT / OUTPUT ================================
// =====================================================

const char* INPUT_FILE = "output.root";
const char* TREE_NAME = "decay_points";
const int MAX_EVENTS = 500000;

// =====================================================
// ===== DATA MODE =====================================
// =====================================================

enum class DataMode {
    Simulation,
    Experimental
};

const DataMode DATA_MODE = DataMode::Simulation;

// =====================================================
// ===== TRACK FINDER ==================================
// =====================================================

enum class TrackingAlgorithm {
    Kalman,
    Hough
};

const TrackingAlgorithm TRACKING_ALGORITHM = TrackingAlgorithm::Kalman;

// =====================================================
// ===== DEBUG =========================================
// =====================================================

const bool DEBUG = false;

// =====================================================
// ===== COMBINATORIAL KALMAN FILTER (CKF) =============
// =====================================================

enum class ZOrder {
    DecreasingZ,
    IncreasingZ
};

const ZOrder CKF_Z_ORDER = ZOrder::DecreasingZ;
const int CKF_MIN_HITS = 3;
const double CKF_CHI2_CUT = 1.0;
const double CKF_PROCESS_NOISE = 0.0;
const double CKF_AMBIGUITY_CHI2 = 4.0;
const bool CKF_GREEDY = true;
const bool CKF_GENEROUS = false;

// =====================================================
// ===== HOUGH TRANSFORM PARAMETERS ====================
// =====================================================

const double HOUGH_A_MIN = -0.09;
const double HOUGH_A_MAX = 0.09;
const double HOUGH_B_MIN = -80.0;
const double HOUGH_B_MAX = 80.0;
const int HOUGH_A_BINS = 6;
const int HOUGH_B_BINS = 13;
const double HOUGH_OVERLAP = 0.03;
const int HOUGH_THRESHOLD = 6;

// =====================================================
// ===== EVENT QUALITY CUTS ============================
// =====================================================

const double CHI2_PER_NDF_CUT = 0.6;
const int EXPECTED_N_TRACKS = 3;

// =====================================================
// ===== OUTPUT ========================================
// =====================================================

const bool PRINT_STATISTICS = true;

#endif // SETTINGS_H
