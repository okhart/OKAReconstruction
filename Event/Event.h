class Event {
public:
    Event() = default;
    ~Event() = default;


    void setHeader(const EventHeader& header);
    const EventHeader& getHeader() const;


    void setTracks(const std::vector<Track>& tracks);
    const std::vector<Track>& getTracks() const;

    void setShowers(const std::vector<Shower>& showers);
    const std::vector<Shower>& getShowers() const;

    void setVertices(const std::vector<Vertex>& vertices);
    const std::vector<Vertex>& getVertices() const;

    void setClusters(const std::vector<Cluster>& clusters);
    const std::vector<Cluster>& getClusters() const;

private:
    EventHeader header;

    std::vector<Track> tracks;
    std::vector<Shower> showers;
    std::vector<Vertex> vertices;
    std::vector<Particle> particles;
};