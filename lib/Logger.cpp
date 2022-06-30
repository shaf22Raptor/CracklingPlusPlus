// Logger Class and Log Buffer Class
#include <Logger.hpp>

using std::string;
using std::ostream;
using std::ofstream;

logBuffer::logBuffer() : isAtStartOfLine(true), fileBuf(NULL), stdBuf(NULL) {};

logBuffer::logBuffer(ostream& stdOut, ostream& destFile) :
    isAtStartOfLine(true),
    fileBuf(NULL),
    stdBuf(NULL)
{
    // Record output stream buffers
    isAtStartOfLine = true;
    fileBuf = destFile.rdbuf();
    stdBuf = stdOut.rdbuf();
}

int logBuffer::sync() {
    // Call sync for both output streams and buffers
    return stdBuf->pubsync() & fileBuf->pubsync();
}

int logBuffer::underflow(int c)
{
    // Do nothing, return EOF
    return EOF;
}

int logBuffer::overflow(int c)
{
    // Insert Timestamp
    if (isAtStartOfLine) {
        time_t rawtime = time(0);
        struct tm* timeinfo = localtime(&rawtime);
        char timestampBuffer[32];
        strftime(timestampBuffer, 32, ">>> %Y-%m-%d %H:%M:%S: ", timeinfo);
        stdBuf->sputn(timestampBuffer, strlen(timestampBuffer));
        fileBuf->sputn(timestampBuffer, strlen(timestampBuffer));
    }
    // Check for new line
    isAtStartOfLine = c == '\n';
    // Add extra newline
    if (isAtStartOfLine) { stdBuf->sputc(c) & fileBuf->sputc(c); }
    // Pass 'c' to output buffers
    return stdBuf->sputc(c) & fileBuf->sputc(c);
}

Logger::Logger(ostream& logSource, string outFile)
{
    // Store Original ref
    origOutputStream = &logSource;
    origSrcBuffer = logSource.rdbuf();

    // Open file stream
    outputFileStream = ofstream(outFile, std::ios::binary);

    // Create custom logBuffer
    customBuffer = logBuffer(*origOutputStream, outputFileStream);

    // Create new ostream
    ostream loggedOutputStream(&customBuffer);

    // Set log source to use new ostream that uses our custom buffer
    origOutputStream->rdbuf(loggedOutputStream.rdbuf());
}

void Logger::close() {
    // Close file stream
    outputFileStream.close();
    // Reset original ref
    origOutputStream->rdbuf(origSrcBuffer);
}