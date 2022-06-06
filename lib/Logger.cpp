// Logger Class and Log Buffer Class
#include <Logger.hpp>

using std::string;
using std::ostream;
using std::ofstream;


logBuffer::logBuffer(ostream& stdOut, ostream& destFile)
{
    isAtStartOfLine = true;
    fileBuf = destFile.rdbuf();
    stdBuf = stdOut.rdbuf();
}

int logBuffer::sync() {
    fileBuf->pubsync();
    return stdBuf->pubsync();
}

int logBuffer::underflow(int c)
{
    return EOF;
}

int logBuffer::overflow(int c)
{
    if (isAtStartOfLine) {
        time_t rawtime = time(0);
        struct tm* timeinfo = localtime(&rawtime);
        char timestampBuffer[32];
        strftime(timestampBuffer, 32, ">>> %Y-%m-%d %H:%M:%S:\t", timeinfo);
        stdBuf->sputn(timestampBuffer, strlen(timestampBuffer));
        fileBuf->sputn(timestampBuffer, strlen(timestampBuffer));
    }
    isAtStartOfLine = c == '\n';

    return stdBuf->sputc(c) & fileBuf->sputc(c);
}

Logger::Logger(ostream& logSource, string outFile)
{
    // Store Original ref
    origOutputStream = &logSource;
    origSrcBuffer = logSource.rdbuf();

    // Open file stream
    outputFileStream = ofstream(outFile);

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