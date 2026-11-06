Event EventBuilder::build() const {
    Event event;

    event.setHeader(header);
    event.setTracks(tracks);
    event.setShowers(showers);
    event.setVertex(vertex);
    event.setParticles(particles);

    return event;
}

void EventBuilder::reset() {
    tracks.clear();
    showers.clear();
    vertex.clear();
    particles.clear();

    header = EventHeader{};
}

void EventBuilder::setTracks(
    const std::vector<Track>& inputTracks)
{
    tracks = inputTracks;
}

void EventBuilder::setShowers(
    const std::vector<Shower>& inputShowers)
{
    showers = inputShowers;
}

void EventBuilder::setVertices(
     Vertex inputVertex)
{
    vertex = inputVertex;
}