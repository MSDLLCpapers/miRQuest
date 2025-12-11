#!/bin/bash

# Docker Runtime Testing Script for miRQuest
# This script tests Docker builds and validates runtime functionality

set -e  # Exit on error

echo "======================================"
echo "miRQuest Docker Build & Runtime Test"
echo "======================================"
echo ""

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Configuration
IMAGE_NAME="mirquest-test"
CONTAINER_NAME="mirquest-test-container"
TEST_PORT=3838
DOCKERFILE="${1:-Dockerfile}"  # Default to Dockerfile, allow override

echo -e "${YELLOW}Configuration:${NC}"
echo "  Dockerfile: $DOCKERFILE"
echo "  Image Name: $IMAGE_NAME"
echo "  Container Name: $CONTAINER_NAME"
echo "  Test Port: $TEST_PORT"
echo ""

# Cleanup function
cleanup() {
    echo ""
    echo "Cleaning up test resources..."
    docker stop "$CONTAINER_NAME" 2>/dev/null || true
    docker rm "$CONTAINER_NAME" 2>/dev/null || true
    if [ "$KEEP_IMAGE" != "true" ]; then
        docker rmi "$IMAGE_NAME" 2>/dev/null || true
    fi
}

# Set trap to cleanup on exit
trap cleanup EXIT

# Step 1: Build the Docker image
echo -e "${YELLOW}Step 1: Building Docker image...${NC}"
echo "Command: docker build -f $DOCKERFILE -t $IMAGE_NAME ."
echo ""

if docker build -f "$DOCKERFILE" -t "$IMAGE_NAME" . ; then
    echo -e "${GREEN}✓ Docker build successful${NC}"
else
    echo -e "${RED}✗ Docker build failed${NC}"
    exit 1
fi

echo ""

# Step 2: Verify image was created
echo -e "${YELLOW}Step 2: Verifying image creation...${NC}"
if docker images | grep -q "$IMAGE_NAME"; then
    IMAGE_SIZE=$(docker images "$IMAGE_NAME" --format "{{.Size}}")
    echo -e "${GREEN}✓ Image created successfully (Size: $IMAGE_SIZE)${NC}"
else
    echo -e "${RED}✗ Image not found${NC}"
    exit 1
fi

echo ""

# Step 3: Start the container
echo -e "${YELLOW}Step 3: Starting container...${NC}"
echo "Command: docker run -d --name $CONTAINER_NAME -p $TEST_PORT:3838 $IMAGE_NAME"

if docker run -d --name "$CONTAINER_NAME" -p "$TEST_PORT:3838" "$IMAGE_NAME"; then
    echo -e "${GREEN}✓ Container started${NC}"
else
    echo -e "${RED}✗ Failed to start container${NC}"
    exit 1
fi

echo ""

# Step 4: Wait for app to initialize
echo -e "${YELLOW}Step 4: Waiting for Shiny app to initialize...${NC}"
echo "This may take 30-60 seconds..."

MAX_WAIT=120  # Maximum wait time in seconds
WAIT_INTERVAL=5
ELAPSED=0

while [ $ELAPSED -lt $MAX_WAIT ]; do
    if docker logs "$CONTAINER_NAME" 2>&1 | grep -q "Listening on"; then
        echo -e "${GREEN}✓ Shiny app is listening${NC}"
        break
    fi

    # Check if container is still running
    if ! docker ps | grep -q "$CONTAINER_NAME"; then
        echo -e "${RED}✗ Container stopped unexpectedly${NC}"
        echo ""
        echo "Container logs:"
        docker logs "$CONTAINER_NAME"
        exit 1
    fi

    echo "  Waiting... ($ELAPSED seconds elapsed)"
    sleep $WAIT_INTERVAL
    ELAPSED=$((ELAPSED + WAIT_INTERVAL))
done

if [ $ELAPSED -ge $MAX_WAIT ]; then
    echo -e "${RED}✗ Timeout waiting for app to start${NC}"
    echo ""
    echo "Container logs:"
    docker logs "$CONTAINER_NAME"
    exit 1
fi

echo ""

# Step 5: Check container health
echo -e "${YELLOW}Step 5: Checking container health...${NC}"
HEALTH_STATUS=$(docker inspect --format='{{.State.Health.Status}}' "$CONTAINER_NAME" 2>/dev/null || echo "none")

if [ "$HEALTH_STATUS" = "none" ]; then
    echo "  No health check configured, skipping..."
elif [ "$HEALTH_STATUS" = "healthy" ]; then
    echo -e "${GREEN}✓ Container is healthy${NC}"
elif [ "$HEALTH_STATUS" = "starting" ]; then
    echo -e "${YELLOW}⚠ Container health is still initializing${NC}"
else
    echo -e "${RED}✗ Container health check failed: $HEALTH_STATUS${NC}"
fi

echo ""

# Step 6: Test HTTP connectivity
echo -e "${YELLOW}Step 6: Testing HTTP connectivity...${NC}"

# Wait a bit more to ensure app is fully ready
sleep 5

if curl -s -o /dev/null -w "%{http_code}" "http://localhost:$TEST_PORT" | grep -q "200"; then
    echo -e "${GREEN}✓ HTTP endpoint responding (200 OK)${NC}"
elif curl -s -o /dev/null -w "%{http_code}" "http://localhost:$TEST_PORT" | grep -q "302"; then
    echo -e "${GREEN}✓ HTTP endpoint responding (302 Redirect)${NC}"
else
    HTTP_CODE=$(curl -s -o /dev/null -w "%{http_code}" "http://localhost:$TEST_PORT")
    echo -e "${RED}✗ HTTP endpoint not responding correctly (HTTP $HTTP_CODE)${NC}"
    echo ""
    echo "Container logs:"
    docker logs "$CONTAINER_NAME" | tail -50
    exit 1
fi

echo ""

# Step 7: Verify key packages are loaded
echo -e "${YELLOW}Step 7: Checking if critical packages loaded...${NC}"
LOGS=$(docker logs "$CONTAINER_NAME" 2>&1)

if echo "$LOGS" | grep -q "miRQuest Docker Runtime"; then
    echo -e "${GREEN}✓ Docker runtime diagnostics executed${NC}"
else
    echo -e "${YELLOW}⚠ Docker runtime diagnostics not found in logs${NC}"
fi

# Check for package loading errors
if echo "$LOGS" | grep -qi "error.*package.*not found"; then
    echo -e "${RED}✗ Package loading errors detected${NC}"
    echo ""
    echo "$LOGS" | grep -i "error.*package"
    exit 1
else
    echo -e "${GREEN}✓ No package loading errors detected${NC}"
fi

echo ""

# Step 8: Display container resource usage
echo -e "${YELLOW}Step 8: Container resource usage:${NC}"
docker stats "$CONTAINER_NAME" --no-stream --format "table {{.Container}}\t{{.CPUPerc}}\t{{.MemUsage}}\t{{.MemPerc}}"

echo ""

# Step 9: Show sample logs
echo -e "${YELLOW}Step 9: Sample application logs:${NC}"
echo "----------------------------------------"
docker logs "$CONTAINER_NAME" 2>&1 | tail -20
echo "----------------------------------------"

echo ""

# Success summary
echo -e "${GREEN}======================================"
echo "✓ All tests passed successfully!"
echo "======================================${NC}"
echo ""
echo "Container Information:"
echo "  Name: $CONTAINER_NAME"
echo "  Port: http://localhost:$TEST_PORT"
echo "  Status: Running"
echo ""
echo "To view live logs:"
echo "  docker logs -f $CONTAINER_NAME"
echo ""
echo "To stop the container:"
echo "  docker stop $CONTAINER_NAME"
echo ""
echo "To keep the image after cleanup, set KEEP_IMAGE=true"
echo ""

# Ask user if they want to keep the container running
read -p "Keep container running for manual testing? (y/N): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    trap - EXIT  # Disable cleanup trap
    echo ""
    echo -e "${GREEN}Container is still running at http://localhost:$TEST_PORT${NC}"
    echo "To stop later: docker stop $CONTAINER_NAME && docker rm $CONTAINER_NAME"
else
    echo ""
    echo "Cleaning up..."
fi
