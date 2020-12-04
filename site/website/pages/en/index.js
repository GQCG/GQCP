/**
 * Copyright (c) Facebook, Inc. and its affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

const React = require('react');

const CompLibrary = require('../../core/CompLibrary.js');

const MarkdownBlock = CompLibrary.MarkdownBlock; /* Used to read markdown */
const Container = CompLibrary.Container;
const GridBlock = CompLibrary.GridBlock;

class HomeSplash extends React.Component {
    render() {
        const { siteConfig, language = '' } = this.props;
        const { baseUrl, docsUrl } = siteConfig;
        const docsPart = `${docsUrl ? `${docsUrl}/` : ''}`;
        const langPart = `${language ? `${language}/` : ''}`;
        const docUrl = (doc) => `${baseUrl}${docsPart}${langPart}${doc}`;

        const SplashContainer = (props) => (
            <div className="homeContainer">
                <div className="homeSplashFade">
                    <div className="wrapper homeWrapper">{props.children}</div>
                </div>
            </div>
        );

        const Logo = (props) => (
            <div className="projectLogo">
                <img src={props.img_src} alt="Project Logo" />
            </div>
        );

        const ProjectTitle = (props) => (
            <h2 className="projectTitle">
                {props.title}
                <small>{props.tagline}</small>
            </h2>
        );

        const PromoSection = (props) => (
            <div className="section promoSection">
                <div className="promoRow">
                    <div className="pluginRowBlock">{props.children}</div>
                </div>
            </div>
        );

        const Button = (props) => (
            <div className="pluginWrapper buttonWrapper">
                <a className="button" href={props.href} target={props.target}>
                    {props.children}
                </a>
            </div>
        );

        return (
            <SplashContainer>
                {/* <Logo img_src={`${baseUrl}img/undraw_monitor.svg`} /> */}
                <div className="inner">
                    <ProjectTitle tagline={siteConfig.tagline} title={siteConfig.title} />
                    <PromoSection>
                        <Button href="https://mybinder.org/v2/gh/GQCG/GQCP/develop?filepath=gqcpy%2Fexamples%2FTry-out-Binder.ipynb">Try It Out!</Button>
                        {/* <Button href={docUrl('doc1.html')}>Example Link</Button>
                        <Button href={docUrl('doc2.html')}>Example Link 2</Button> */}
                    </PromoSection>
                </div>
            </SplashContainer>
        );
    }
}

class Index extends React.Component {
    render() {
        const { config: siteConfig, language = '' } = this.props;
        const { baseUrl } = siteConfig;

        const Block = (props) => (
            <Container
                padding={['bottom', 'top']}
                id={props.id}
                background={props.background}>
                <GridBlock
                    align={props.align}
                    contents={props.children}
                    layout={props.layout}
                />
            </Container>
        );

        const SustainableResearchTest = () => (
            <div
                background="dark"
                className="productShowcaseSection paddingBottom"
                style={{ textAlign: 'left' }}>
                <h2>Sustainable development & reproducible research</h2>
                <MarkdownBlock>
                    At GQCG, our activities are centered around electronic structure for molecules. We're trying to help create a community of quantum chemists that have a common mindset, suitable for (academic) research in the 21st century. We're focused on sustainable and reproducible research, which is why our software developments are all open-source.

                    In order to support our research and development, GQCP saw its first light in 2017 and has been growing ever since.
                </MarkdownBlock>
            </div>
        );

        const SustainableResearch = () => (
            <Block background="dark">
                {[
                    {
                        content:
                            'At the Ghent Quantum Chemistry Group, GQCG for short, our activities are centered around electronic structure for molecules. We are trying to help create a community of quantum chemists that have a common mindset, suitable for (academic) research in the 21st century. We are focused on sustainable and reproducible research, which is why our software developments are all open-source.' +
                            'In order to support our research and development, GQCP saw its first light in 2017 and has been growing ever since.',
                        image: `${baseUrl}img/gqcg_logo-01.png`,
                        imageAlign: 'right',
                        title: 'Sustainable development & reproducible research',
                    },
                ]}
            </Block>
        );

        // const Description = () => (
        //     <Block background="light">
        //         {[
        //             {
        //                 content:
        //                     'Knowdes, a portmanteau of knowledge and nodes, is a collection of interconnected knowledge nodes. It is a web of knowledge that is aimed at supporting the research performed at GQCG.' +
        //                     'At GQCG, Knowdes is one of our premier sources of information. Everything from research notes to theoretical support and elaborations can be found there.' +
        //                     'You can find Knowdes [here](https://gqcg-res.github.io/knowdes/index.html).',
        //                 image: `${baseUrl}img/knowdes.jpg`,
        //                 imageAlign: 'left|middle',
        //                 title: 'Knowdes',
        //             },
        //         ]}
        //     </Block>
        // );

        const Blocks = () => (
            <Block background="light">
                {[
                    {
                        content:
                            'Conceptually, GQCP is like Lego. Much like our favorite toys, GQCP provides building blocks to experiment with. In just a handful lines of code, users can perform complex computations and developers can synthesize new and modify existing electronic structure methods. Curious? Try it out!',
                        image: `${baseUrl}img/blocks.png`,
                        imageAlign: 'left',
                        title: 'GQCWhy?',
                    },
                ]}
            </Block>
        );

        const Features = () => (
            <Block
                layout="twoColumn"
                align='center'>
                {[
                    {
                        title: "Python bindings",
                        image: `${baseUrl}img/bindings.png`,
                        imageAlign: 'top',
                        content: "We use [pybind11](https://github.com/pybind/pybind11) to generate Python bindings for our C++ library. Using `gqcpy` as a Python module, we can embrace Python's present role as a data manipulating language. Gone are the days of providing input files or writing executables, with GQCPY and Jupyter Notebooks you can naturally adapt a work flow that is both playful and systematic at the same time."
                    },
                    {
                        title: 'C++ library',
                        image: `${baseUrl}img/python-and-cpp.png`,
                        imageAlign: 'top',
                        content: "GQCP is modern at its core. It is natively written in C++, so we have access to state-of-the-art software techniques and compilers."
                    },
                ]}
            </Block>
        );

        const Showcase = () => {
            if ((siteConfig.references || []).length === 0) {
                return null;
            }

            const showcase = siteConfig.references
                .filter((reference) => reference.pinned)
                .map((reference) => (
                    <a href={reference.infoLink} key={reference.infoLink}>
                        {reference.caption}
                    </a>
                ));

            const pageUrl = (page) =>
                baseUrl + (language ? `${language}/` : '') + page;

            return (
                <div className="productShowcaseSection paddingBottom">
                    <h2>Academic references</h2>
                    <p>This publication, amongst others, used GQCP.</p>
                    <div className="logos">{showcase}</div>
                    <div className="more-users">
                        <a className="button" href={pageUrl('references.html')}>
                            All publications that used {siteConfig.title}
                        </a>
                    </div>
                </div>
            );
        };

        const Developers = () => {
            if ((siteConfig.developers || []).length === 0) {
                return null;
            }

            const showcase = siteConfig.developers.map((developer) => (
                <a href={developer.infoLink} key={developer.infoLink}>
                    <img src={developer.image} alt={developer.caption} title={developer.caption} />
                </a>
            ));

            return (
                <div className="productShowcaseSection paddingBottom">
                    <h2>The GQCP development team</h2>
                    <div className="logos">{showcase}</div>
                </div>
            );
        };


        return (
            <div>
                <HomeSplash siteConfig={siteConfig} language={language} />
                <div className="mainContainer">
                    <Features />
                    <Blocks />
                    <SustainableResearch />
                    <Showcase />
                    <Developers />
                </div>
            </div>
        );
    }
}

module.exports = Index;
