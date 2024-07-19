const midi = require("midi");
const sharp = require("sharp");
const fs = require("fs")

const width = 160;
const height = 80;

function intToTwoByteHex(int) {
  // Ensure the integer is within the range of two bytes (0 to 65535)
  if (int < 0 || int > 65535) {
    throw new RangeError("Input integer must be between 0 and 65535.");
  }

  // Convert the integer to a hexadecimal string and pad it to 4 characters
  let hexString = int.toString(16).padStart(4, '0');

  // Split the hex string into two bytes
  let highByte = parseInt(hexString.slice(0, 2), 16);
  let lowByte = parseInt(hexString.slice(2, 4), 16);

  return [highByte, lowByte];
}

function magicNumber(seed) {
  function magicFormula(n) { return 128 * n - 8 }

  let n = 1
  while (magicFormula(n) < seed) n++

  return magicFormula(n - 1)
}

function encodePng(originalBuffer) {
  const convertedBuffer = [];

  for (let i = 0; i < originalBuffer.length; i += 7) {
    const group = originalBuffer.slice(i, i + 7);
    let controlByte = 0;

    for (let j = 0; j < group.length; j++) {
      if (group[j] >= 128) {
        controlByte |= 1 << j;
        group[j] -= 128;
      }
    }

    convertedBuffer.push(controlByte, ...group);
  }

  return Buffer.from(convertedBuffer);
}

function generateMessageMetadata(x, y, decodedLength, encodedLength) {
  const imageCoordinates = [x, 0x00, y, 0x00];

  const messageLength = encodedLength + 16

  const magic = magicNumber(encodedLength);
  const meta1Dec = messageLength + magic
  const meta1 = intToTwoByteHex(meta1Dec)

  let meta2 = [0x00, 0x20]
  let meta3 = intToTwoByteHex(decodedLength).toReversed()

  if (meta3[0] >= 128) {
    meta2 = [0x20, 0x20]
    meta3 = intToTwoByteHex(decodedLength - 128).toReversed()
  }

  return [
    ...meta1,
    ...meta2,
    ...imageCoordinates,
    ...meta3,
  ];
}

function generateMessage(imageData, x, y) {
  const messageHeader = [0xf0, 0x47, 0x7f, 0x4a, 0x04];
  const encoded = encodePng(imageData)
  const messageMetadata = generateMessageMetadata(x, y, imageData.length, encoded.length)
  const terminator = [0xf7]

  return [
    ...messageHeader,
    ...messageMetadata,
    ...encoded,
    ...terminator
  ]
}

function sendChunks(output, chunks) {
  output.sendMessage(generateMessage(chunks[0], 0, 0))
  output.sendMessage(generateMessage(chunks[1], 0, 60))
  output.sendMessage(generateMessage(chunks[2], 60, 0))
  output.sendMessage(generateMessage(chunks[3], 60, 60))
  output.sendMessage(generateMessage(chunks[4], 120, 0))
  output.sendMessage(generateMessage(chunks[5], 120, 60))
}

async function generateImage() {
  const pixels = Buffer.alloc(width * height * 3);
  for (let i = 0; i < pixels.length; i++) {
    pixels[i] = Math.floor(Math.random() * 256);
  }
  const image = sharp(pixels, {
    raw: {
      width: width,
      height: height,
      channels: 3
    }
  })
  return await image.toFormat("png").toBuffer();
}

function splitImage(image) {
  const splitCoordinates = [
    [0, 0, 60, 60],
    [0, 60, 60, 20],
    [60, 0, 60, 60],
    [60, 60, 60, 20],
    [120, 0, 40, 60],
    [120, 60, 40, 20],
  ];

  return Promise.all(
    splitCoordinates.map(async (splitCoordinate, i) => {
      const [left, top, width, height] = splitCoordinate;
      const splitImage = sharp(image).extract({ left, top, width, height });
      return await splitImage
        .toFormat("png", {
          compressionLevel: 9,
          colours: 8,
        }).toBuffer()
    })
  );
}

(async () => {
  const output = new midi.Output()

  let port = 0
  for (let i = 0; i < output.getPortCount(); i++) {
    let portName = output.getPortName(i)
    if (portName.includes('MPC Studio')) {
      port = i
      break
    }
  }
  console.log(port, output.getPortName(port))
  output.openPort(port)

  setInterval(async () => {
    const image = await generateImage()
    const imageChunks = await splitImage(image)
    sendChunks(output, imageChunks)
  }, 100)

  // output.closePort()
})()
