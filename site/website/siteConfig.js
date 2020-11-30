/**
 * Copyright (c) Facebook, Inc. and its affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

// See https://docusaurus.io/docs/site-config for all the possible
// site configuration options.

// List of projects/orgs using your project for the users page.
const references = [
    {
        caption: 'Hierarchies of quantum chemical descriptors induced by statistical analyses of domain occupation number operators',
        // You will need to prepend the image path with your baseUrl
        // if it is not '/', like: '/test-site/img/image.jpg'.
        // image: '/img/undraw_open_source.svg',
        infoLink: 'https://onlinelibrary.wiley.com/doi/abs/10.1002/wcms.1456',
        pinned: true,
    },
];

const developers = [
    {
        caption: '@guacke',
        // You will need to prepend the image path with your baseUrl
        // if it is not '/', like: '/test-site/img/image.jpg'.
        image: '/GQCP/img/guacke.png',
        infoLink: 'https://github.com/guacke',
    },

    {
        caption: '@lelemmen',
        // You will need to prepend the image path with your baseUrl
        // if it is not '/', like: '/test-site/img/image.jpg'.
        image: '/GQCP/img/lelemmen.jpeg',
        infoLink: 'https://github.com/lelemmen',
    },

    {
        caption: '@dariatols',
        // You will need to prepend the image path with your baseUrl
        // if it is not '/', like: '/test-site/img/image.jpg'.
        image: '/GQCP/img/dariatols.png',
        infoLink: 'https://github.com/dariatols',
    },

    {
        caption: '@xdvriend',
        // You will need to prepend the image path with your baseUrl
        // if it is not '/', like: '/test-site/img/image.jpg'.
        image: '/GQCP/img/xdvriend.JPG',
        infoLink: 'https://github.com/xdvriend',
    },
];

const siteConfig = {
    title: 'GQCP', // Title for your website.
    tagline: 'The Ghent Quantum Chemistry Package',
    url: 'https://GQCG.github.io', // Your website URL
    baseUrl: '/GQCP/', // Base URL for your project */
    // For github.io type URLs, you would set the url and baseUrl like:
    //   url: 'https://facebook.github.io',
    //   baseUrl: '/test-site/',

    // Used for publishing and more
    projectName: 'GQCP',
    organizationName: 'GQCG',
    // For top-level user or org sites, the organization is still the same.
    // e.g., for the https://JoelMarcey.github.io site, it would be set like...
    //   organizationName: 'JoelMarcey'

    // For no header links in the top nav bar -> headerLinks: [],
    headerLinks: [
        { doc: 'documentation', label: 'Documentation' },
        { doc: 'api', label: 'API' },
        { page: 'help', label: 'Help' },
        { blog: true, label: 'Blog' },
    ],

    // If you have references set above, you add it here:
    references,

    // If you have developers set above, you add it here:
    developers,

    /* Path to images for header/footer. Uncomment if you want to enable it. */
    // headerIcon: 'img/gqcg_logo-01.png',
    // footerIcon: 'img/gqcg_logo-01.png',
    // favicon: 'img/gqcg_logo-01.png',

    /* Colors for website */
    colors: {
        primaryColor: '#1E64C8',
        secondaryColor: '#528ad9',
    },

    /* Custom fonts for website */
    /*
    fonts: {
      myFont: [
        "Times New Roman",
        "Serif"
      ],
      myOtherFont: [
        "-apple-system",
        "system-ui"
      ]
    },
    */

    // This copyright info is used in /core/Footer.js and blog RSS/Atom feeds.
    copyright: `Copyright © ${new Date().getFullYear()} GQCG`,

    highlight: {
        // Highlight.js theme to use for syntax highlighting in code blocks.
        theme: 'default',
    },

    // Add custom scripts here that would be placed in <script> tags.
    scripts: ['https://buttons.github.io/buttons.js'],

    // On page navigation for the current documentation page.
    onPageNav: 'separate',
    // No .html extensions for paths.
    cleanUrl: true,

    // Open Graph and Twitter card images.
    ogImage: 'img/undraw_online.svg',
    twitterImage: 'img/undraw_tweetstorm.svg',

    // For sites with a sizable amount of content, set collapsible to true.
    // Expand/collapse the links and subcategories under categories.
    // docsSideNavCollapsible: true,

    // Show documentation's last contributor's name.
    // enableUpdateBy: true,

    // Show documentation's last update time.
    // enableUpdateTime: true,

    // You may provide arbitrary config keys to be used as needed by your
    // template. For example, if you need your repo's URL...
    // repoUrl: 'https://github.com/facebook/test-site',
};

module.exports = siteConfig;
