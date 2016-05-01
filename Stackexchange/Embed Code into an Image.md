## Embedding Code into an Image ##

*(Thank god there is Mr. McLoone with his creative mind)*

Say for instance we have a plot:

    Plot[Sin[x], {x, 0, 2 Pi}]

![enter image description here][1]


To generate an image file:

    Export["~/plot.png", Plot[Sin[x], {x, 0, 2 Pi}]];
    image = Import["~/plot.png"]; DeleteFile["~/plot.png"];

(is there a shortcut for this?)

**Edit** (thank you to rm -rf)

miss the wood for the trees:

    image = Rasterize[Plot[Sin[x], {x, 0, 2 Pi}]];

**Edit end**

If we say now, that we use for every color channel the least significant bit to embed the code, we first have to create a truncated variant of the original
image forcing strictly 8-bits per channel:

    truncImage = BitAnd[ImageData[image, "Byte"], 2^^11111110];

Next we have to convert the code into a sequence of bits and insert each of them into those empty bits; padding the rest with zeros:

    code = PadRight[
        Flatten[IntegerDigits[ToCharacterCode@Compress["Plot[Sin[x],{x,0,2Pi}]"], 2,8]],
        Apply[Times, Dimensions[ImageData[image, "Byte"]]]];

Now we have to merge the truncated image with *code*:

    codeImage = Image[truncImage + Fold[Partition, code,
     Reverse@Rest[Dimensions[ImageData[image, "Byte"]]]], "Byte"];

Now the code is embedded into that image. In order to extract the code, we've to reverse the process:

    secCode = FromDigits[#, 2] & /@ (Partition[
     Flatten@BitAnd[ImageData[codeImage, "Byte"], 1], 8]);

and then:

    Uncompress@FromCharacterCode[secCode]

=> Plot[Sin[x],{x,0,2Pi}]

or:

    Uncompress@FromCharacterCode[secCode] // ToExpression

![enter image description here][1]

 [1]: http://i.stack.imgur.com/HdRd5.png

 
