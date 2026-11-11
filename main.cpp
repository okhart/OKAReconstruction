int main() {

    EventBuilder eventBuilder;

    while (/* there are events */) {

        eventBuilder.reset();

        // Получение результатов реконструкции
        auto tracks = trackReconstruction.run();
        auto showers = showerReconstruction.run();
        auto vertices = vertexReconstruction.run();
        auto particles = particleIdentification.run();
        auto clusters = clusterReconstruction.run();

        // Передача результатов EventBuilder
        eventBuilder.setTracks(tracks);
        eventBuilder.setShowers(showers);
        eventBuilder.setVertices(vertices);
        eventBuilder.setParticles(particles);
        eventBuilder.setClusters(clusters);

        // Формирование итогового события
        Event event = eventBuilder.build();

        // Передача события следующему модулю
    }

    return 0;
}