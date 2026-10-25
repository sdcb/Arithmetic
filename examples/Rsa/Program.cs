// Textbook RSA built from BigInt primitives: key generation (NextPrime, ModInverse),
// encryption/decryption (ModPow, with the CRT speed-up) and PKCS#1 v1.5 signatures that
// interoperate with System.Security.Cryptography.RSA in both directions.
// FOR LEARNING ONLY: no constant-time arithmetic, no padding for encryption.
// Usage: dotnet run --project examples/Rsa -c Release -- [modulus-bits]
using System.Diagnostics;
using System.Security.Cryptography;
using System.Text;
using Sdcb.Arithmetic;
using Sdcb.Arithmetic.Examples;

PerfWarning.Print();

int bits = args.Length > 0 ? int.Parse(args[0]) : 2048;
var sw = Stopwatch.StartNew();

// --- Key generation ---------------------------------------------------------------
BigInt e = 65537;
BigInt p, q, n, d;
do
{
    p = RandomPrime(bits / 2);
    q = RandomPrime(bits - bits / 2);
    n = p * q;
} while (p == q || n.GetBitLength() != bits || BigInt.Gcd(e, (p - 1) * (q - 1)) != 1);
d = BigInt.ModInverse(e, BigInt.Lcm(p - 1, q - 1));
BigInt dp = d % (p - 1), dq = d % (q - 1), qInv = BigInt.ModInverse(q, p);
Console.WriteLine($"generated {bits}-bit key in {sw.ElapsedMilliseconds} ms");
Console.WriteLine($"n = {n.ToString(16)[..32]}... ({n.GetBitLength()} bits)");

// --- Encrypt / decrypt ---------------------------------------------------------------
string message = "Hello, Sdcb.Arithmetic!";
BigInt m = BigInt.FromMagnitudeBytes(Encoding.UTF8.GetBytes(message), bigEndian: true);
BigInt c = BigInt.ModPow(m, e, n);
BigInt plain = BigInt.ModPow(c, d, n), crt = DecryptCrt(c); // warm up
sw.Restart();
for (int i = 0; i < 20; i++) plain = BigInt.ModPow(c, d, n);
double plainMs = sw.Elapsed.TotalMilliseconds / 20;
sw.Restart();
for (int i = 0; i < 20; i++) crt = DecryptCrt(c);
double crtMs = sw.Elapsed.TotalMilliseconds / 20;
Console.WriteLine($"decrypted: \"{Encoding.UTF8.GetString(plain.ToMagnitudeBytes(bigEndian: true))}\"");
Console.WriteLine($"decrypt: c^d mod n {plainMs:F2} ms, CRT {crtMs:F2} ms (same result: {crt == plain})");

// --- Signatures, cross-checked with .NET's RSA ---------------------------------------
using RSA rsa = RSA.Create();
rsa.ImportParameters(ToRsaParameters());
byte[] data = Encoding.UTF8.GetBytes("sign me");

byte[] ours = Sign(data);
bool dotnetVerifiesOurs = rsa.VerifyData(data, ours, HashAlgorithmName.SHA256, RSASignaturePadding.Pkcs1);
byte[] theirs = rsa.SignData(data, HashAlgorithmName.SHA256, RSASignaturePadding.Pkcs1);
bool weVerifyTheirs = Verify(data, theirs);
Console.WriteLine($".NET verifies our signature: {dotnetVerifiesOurs}");
Console.WriteLine($"we verify .NET's signature : {weVerifyTheirs}");
Console.WriteLine($"signatures identical       : {ours.AsSpan().SequenceEqual(theirs)}");
return dotnetVerifiesOurs && weVerifyTheirs ? 0 : 1;

// A random prime with exactly `bits` bits (top two bits set so p*q has the full length).
static BigInt RandomPrime(int bits)
{
    byte[] buf = RandomNumberGenerator.GetBytes((bits + 7) / 8);
    BigInt x = BigInt.FromMagnitudeBytes(buf) >> (buf.Length * 8 - bits);
    x |= BigInt.One << (bits - 1) | BigInt.One << (bits - 2);
    return x.NextPrime();
}

BigInt DecryptCrt(BigInt c)
{
    BigInt m1 = BigInt.ModPow(c, dp, p), m2 = BigInt.ModPow(c, dq, q);
    BigInt h = qInv * (m1 - m2) % p;
    if (h.Sign < 0) h += p;
    return m2 + h * q;
}

// EMSA-PKCS1-v1_5 with SHA-256: 00 01 FF..FF 00 || DigestInfo || hash.
BigInt EncodePkcs1(byte[] data)
{
    byte[] prefix = Convert.FromHexString("3031300d060960864801650304020105000420");
    byte[] em = new byte[(bits + 7) / 8];
    byte[] hash = SHA256.HashData(data);
    em[1] = 0x01;
    int tail = prefix.Length + hash.Length;
    em.AsSpan(2, em.Length - 3 - tail).Fill(0xFF);
    prefix.CopyTo(em, em.Length - tail);
    hash.CopyTo(em, em.Length - hash.Length);
    return BigInt.FromMagnitudeBytes(em, bigEndian: true);
}

byte[] Sign(byte[] data) => Fixed(DecryptCrt(EncodePkcs1(data)), (bits + 7) / 8);

bool Verify(byte[] data, byte[] signature) =>
    BigInt.ModPow(BigInt.FromMagnitudeBytes(signature, bigEndian: true), e, n) == EncodePkcs1(data);

RSAParameters ToRsaParameters()
{
    int len = (bits + 7) / 8, half = (len + 1) / 2;
    return new RSAParameters
    {
        Modulus = Fixed(n, len), Exponent = e.ToMagnitudeBytes(bigEndian: true), D = Fixed(d, len),
        P = Fixed(p, half), Q = Fixed(q, half), DP = Fixed(dp, half), DQ = Fixed(dq, half), InverseQ = Fixed(qInv, half),
    };
}

// Big-endian magnitude left-padded to a fixed length.
static byte[] Fixed(BigInt v, int length)
{
    byte[] raw = v.ToMagnitudeBytes(bigEndian: true);
    byte[] result = new byte[length];
    raw.CopyTo(result, length - raw.Length);
    return result;
}
