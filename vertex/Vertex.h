#include "../tracking/track.h"
#include "../common/Point.h"
class Vertex : public Point{
    
};

enum class VertexStatus : uint8_t {
    Valid       = 0,
    NotConverged = 1,   // фиттер не сошёлся
    TooFewTracks = 2,   // недостаточно треков
    Chi2TooHigh  = 3,   // Chi2 выше порога
    Rejected     = 4    // отброшен
};


struct RecoVertex {
    Point                   position;        
    double                  chi2      = 0.0; // chi2 фиттера
    VertexType              type     = VertexType::Unknown;
    VertexStatus            status   = VertexStatus::Valid;

    // Индексы треков
    std::vector<track>   tracks;


    bool valid() const noexcept {
        return status == VertexStatus::Valid;
    }
};

// Кандидат на вершину
class VertexCandidate {
public:
    VertexCandidate();
private:
    Point                   pos;
    std::vector<track>   tracks;
};

// Абстрактный алгоритм поиска вершины
class VertexFinder {
public:
    virtual ~VertexFinder() = default;

    virtual std::vector<VertexCandidate> find(
        std::vector<track>              tracks,
        const std::optional<Point>&   primary_vertex = std::nullopt) = 0;

    virtual const char* name() const noexcept = 0;
};

