# Stage 1: The Builder
# This stage installs build tools and compiles the C++ code.
FROM python:3.10.12 as builder

WORKDIR /build

# Install C++ build dependencies (add libsqlite3-dev)
RUN apt-get update && \
    apt-get install -y --no-install-recommends g++ libeigen3-dev nlohmann-json3-dev libsqlite3-dev && \
    rm -rf /var/lib/apt/lists/*

# Copy only the C++ source code needed for compilation
COPY src/cpp/ /build/src/cpp/

# Compile the C++ shared library (link against sqlite3)
RUN g++ -I/usr/include/eigen3 -I/usr/include/nlohmann -shared -o /build/interface.so /build/src/cpp/interface.cpp /build/src/cpp/FilterStack.cpp /build/src/cpp/core.cpp /build/src/cpp/input.cpp -fPIC -fopenmp -lsqlite3


# Stage 2: The Final Image
# This is the lean image that will actually run.
FROM python:3.10.12

WORKDIR /app

# Install runtime dependencies (no build tools needed here)
RUN apt-get update && apt-get install -y --no-install-recommends libsqlite3-dev && rm -rf /var/lib/apt/lists/*

# Copy Python requirements and install them first to leverage Docker cache
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
RUN pip install rq

# Copy the rest of your application code
COPY . /app

# Copy the compiled .so file from the builder stage, overwriting the one from the previous copy
COPY --from=builder /build/interface.so /app/src/cpp/interface.so

# Set Environment Variables
ENV FLASK_APP=gui
ENV PYTHONPATH=/app
ENV FLASK_ENV=development

# Expose port and run the application
EXPOSE 5000
CMD ["flask", "run", "--host=0.0.0.0"]
