/**
 * Copyright (c) Facebook, Inc. and its affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

const React = require('react');

class Footer extends React.Component {
    docUrl(doc) {
        const baseUrl = this.props.config.baseUrl;
        const docsUrl = this.props.config.docsUrl;
        const docsPart = `${docsUrl ? `${docsUrl}/` : ''}`;
        return `${baseUrl}${docsPart}${doc}`;
    }

    render() {
        return (
            <footer className="nav-footer" id="footer">
                <section className="sitemap">
                    <a href={this.props.config.baseUrl} className="nav-home">
                        {this.props.config.footerIcon && (
                            <img
                                src={this.props.config.baseUrl + this.props.config.footerIcon}
                                alt={this.props.config.title}
                                width="66"
                                height="58"
                            />
                        )}
                    </a>
                    <div>
                        <h5>Docs</h5>
                        <a href={this.docUrl('user_documentation/user_getting_started.html')}>
                            User documentation
            </a>
                        <a href={this.docUrl('developer_documentation/developer_getting_started.html')}>Developer documentation</a>
                        <a href={this.docUrl('doc4.html')}>
                            API Reference
            </a>
                    </div>
                    <div>
                        <h5>Community</h5>
                        <a href={`${this.props.config.baseUrl}users`}>Academic references</a>
                        <a
                            href="https://github.com/GQCG/GQCP/graphs/contributors"
                            target="_blank"
                            rel="noreferrer noopener">
                            Contributors
            </a>

                    </div>
                    <div>
                        <h5>More</h5>
                        <a href={`${this.props.config.baseUrl}blog`}>Blog</a>
                        <a href="https://github.com/GQCG/GQCP">GitHub</a>
                        <a
                            className="github-button"
                            href="https://github.com/GQCG/GQCP"
                            data-icon="octicon-star"
                            data-count-href="/stargazers"
                            data-show-count="true"
                            data-count-aria-label="# stargazers on GitHub"
                            aria-label="Star this project on GitHub">
                            Star GQCP
            </a>

                    </div>
                </section>

                <section className="copyright">{this.props.config.copyright}</section>
            </footer>
        );
    }
}

module.exports = Footer;
